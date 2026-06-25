#!/usr/bin/env python3
"""
critique_cli.py — standalone LLM code review for Hermes workflows.

Usage:
  python3 scripts/critique_cli.py pool_seq rnaseq   # review specific
  python3 scripts/critique_cli.py                   # review all 11
  python3 scripts/critique_cli.py --changed         # only git-modified workflows
"""

import argparse
import asyncio
import os
import subprocess
import sys
from pathlib import Path

try:
    import anthropic
except ImportError:
    sys.exit("ERROR: anthropic SDK not installed. Run: pip install anthropic")

HERMES_DIR    = Path(__file__).resolve().parent.parent
WF_DIR        = HERMES_DIR / "workflows"
ALL_WORKFLOWS = [
    "pool_seq", "snp_association", "genomic_prediction",
    "atacseq", "scrnaseq", "metagenome", "spatial", "report",
    "rnaseq", "wgs_snp", "genome_annotation",
]

ISSUE_SCHEMA = {
    "type": "object",
    "properties": {
        "workflow": {"type": "string"},
        "critical_issues": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "category":      {"type": "string"},
                    "location":      {"type": "string"},
                    "description":   {"type": "string"},
                    "suggested_fix": {"type": "string"},
                },
                "required": ["category", "location", "description", "suggested_fix"],
            },
        },
        "minor_issues": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "category":      {"type": "string"},
                    "location":      {"type": "string"},
                    "description":   {"type": "string"},
                    "suggested_fix": {"type": "string"},
                },
                "required": ["category", "location", "description", "suggested_fix"],
            },
        },
        "verdict": {
            "type": "string",
            "enum": ["ready_to_run", "minor_fixes_needed", "critical_fixes_required"],
        },
        "summary": {"type": "string"},
    },
    "required": ["workflow", "critical_issues", "minor_issues", "verdict", "summary"],
}


# ── helpers ───────────────────────────────────────────────────────────────────

def get_changed_workflows() -> list[str]:
    """Return workflow names that changed in the last git commit."""
    try:
        out = subprocess.check_output(
            ["git", "diff", "HEAD~1", "--name-only"],
            cwd=str(HERMES_DIR), text=True, stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        return []
    names: set[str] = set()
    for line in out.splitlines():
        if line.startswith("workflows/"):
            parts = line.split("/")
            if len(parts) >= 2:
                names.add(parts[1])
    return [n for n in names if (WF_DIR / n).exists()]


def bundle_files(wf_dir: Path) -> str:
    parts = []
    for target, glob_pat in [
        (wf_dir / "Snakefile",            None),
        (wf_dir / "config_template.yaml", None),
        (wf_dir / "envs",                 "*.yaml"),
        (wf_dir / "scripts",              "*.R"),
        (wf_dir / "scripts",              "*.py"),
    ]:
        if glob_pat:
            if target.exists():
                for f in sorted(target.glob(glob_pat)):
                    parts.append(f"=== {f.relative_to(wf_dir)} ===\n{f.read_text(errors='replace')}")
        elif target.exists():
            parts.append(f"=== {target.name} ===\n{target.read_text(errors='replace')}")
    return "\n\n".join(parts)


def build_prompt(name: str, files_content: str, checklist: str) -> str:
    return f"""You are a senior bioinformatics engineer doing an independent code review of the
"{name}" Snakemake workflow. Review with fresh eyes. Find bugs that cause crashes or wrong
results in production. Cost of a missed bug: days of lost compute on a 160-core server.

## Workflow files
{files_content}

## PIPELINE_CHECKLIST.md
{checklist}

## What to check
**Docker (CRITICAL):** --user $(id -u):$(id -g) on every docker run (EXCEPTION: BRAKER3 root).
mkdir -p output BEFORE docker run. busybox chown for root-owned outputs.

**Shell blocks (CRITICAL):** set -euo pipefail first line. test -s {{output}} after every command.
sed back-references: use \\\\1 not \\1 (Python processes Snakemake shell strings first).

**samtools (CRITICAL):** fixmate requires name-sorted input (sort -n before fixmate).

**Index files:** .fai/.bai/.csi/.tbi declared as inputs where tools require them.
samtools >= 1.12 writes .csi for non-standard references.

**Conda envs:** every library()/import must have a matching package in the conda yaml.

**Statistics:** appropriate normalization, multiple-testing correction, empty-output guards.

Use the submit_review tool. CRITICAL = crash or wrong result. MINOR = style/sub-optimal.
Include file:line for every issue."""


# ── API calls ─────────────────────────────────────────────────────────────────

async def critique_one(client: anthropic.AsyncAnthropic, name: str, checklist: str) -> dict | None:
    wf_dir = WF_DIR / name
    if not wf_dir.exists():
        print(f"[{name}] SKIP: directory not found", flush=True)
        return None
    files_content = bundle_files(wf_dir)
    print(f"[{name}] reviewing ({len(files_content)//1000}k chars) …", flush=True)
    try:
        response = await client.messages.create(
            model       = "claude-sonnet-4-6",
            max_tokens  = 4096,
            tools       = [{"name": "submit_review",
                            "description": "Submit structured code review findings",
                            "input_schema": ISSUE_SCHEMA}],
            tool_choice = {"type": "tool", "name": "submit_review"},
            messages    = [{"role": "user",
                            "content": build_prompt(name, files_content, checklist)}],
        )
    except anthropic.APIError as exc:
        print(f"[{name}] API error: {exc}", flush=True)
        return None
    for block in response.content:
        if block.type == "tool_use" and block.name == "submit_review":
            result = dict(block.input)
            result["workflow"] = name
            return result
    return None


async def synthesize(client: anthropic.AsyncAnthropic, results: list[dict]) -> str:
    details = []
    for r in results:
        crit = "\n".join(
            f"  [{i['category']}] {i['location']}\n    {i['description']}\n    FIX: {i['suggested_fix']}"
            for i in r.get("critical_issues", [])
        ) or "  (none)"
        minor = "\n".join(
            f"  [{i['category']}] {i['location']}: {i['description']}"
            for i in r.get("minor_issues", [])
        ) or "  (none)"
        details.append(
            f"=== {r['workflow']} [{r['verdict']}] ===\n"
            f"CRITICAL ({len(r.get('critical_issues', []))}):\n{crit}\n"
            f"MINOR ({len(r.get('minor_issues', []))}):\n{minor}\n"
            f"Summary: {r.get('summary', '')}"
        )
    response = await client.messages.create(
        model      = "claude-sonnet-4-6",
        max_tokens = 2048,
        messages   = [{"role": "user", "content": (
            "Synthesize these code-review results for Hermes bioinformatics workflows.\n\n"
            + "\n\n".join(details)
            + "\n\n"
            "## 1. Overall verdict\n## 2. Workflows ready to run\n"
            "## 3. Critical fixes (by workflow; exact code)\n"
            "## 4. Systemic patterns\n## 5. Next steps\n"
            "Be specific. No filler."
        )}],
    )
    return response.content[0].text


# ── main ──────────────────────────────────────────────────────────────────────

async def main(targets: list[str]) -> int:
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        sys.exit("ERROR: ANTHROPIC_API_KEY not set")

    checklist = (HERMES_DIR / "PIPELINE_CHECKLIST.md").read_text() \
        if (HERMES_DIR / "PIPELINE_CHECKLIST.md").exists() else ""

    print(f"\n{'='*60}", flush=True)
    print(f"Hermes Critique — reviewing: {', '.join(targets)}", flush=True)
    print(f"{'='*60}\n", flush=True)

    client  = anthropic.AsyncAnthropic(api_key=api_key)
    raw     = await asyncio.gather(*[critique_one(client, n, checklist) for n in targets],
                                   return_exceptions=True)
    results = [r for r in raw if isinstance(r, dict)]

    # ── per-workflow summary ──────────────────────────────────────────────────
    print(f"\n{'─'*60}", flush=True)
    has_critical = False
    for r in results:
        nc = len(r.get("critical_issues", []))
        nm = len(r.get("minor_issues",    []))
        verdict = r["verdict"]
        if verdict == "critical_fixes_required":
            has_critical = True
        icon = "✗" if verdict == "critical_fixes_required" else ("!" if verdict == "minor_fixes_needed" else "✓")
        print(f"  {icon}  {r['workflow']:25s}  {verdict}  ({nc} critical, {nm} minor)", flush=True)

        for issue in r.get("critical_issues", []):
            print(f"        CRITICAL [{issue['category']}] {issue['location']}", flush=True)
            print(f"          {issue['description']}", flush=True)
            print(f"          FIX: {issue['suggested_fix']}", flush=True)

    n_critical = sum(len(r.get("critical_issues", [])) for r in results)
    n_minor    = sum(len(r.get("minor_issues",    [])) for r in results)
    print(f"\nTotal: {len(results)} reviewed  |  {n_critical} critical  |  {n_minor} minor", flush=True)

    # ── synthesis ─────────────────────────────────────────────────────────────
    if results:
        print(f"\n{'='*60}", flush=True)
        print("SYNTHESIS", flush=True)
        print(f"{'='*60}", flush=True)
        syn = await synthesize(client, results)
        print(syn, flush=True)

    print(f"\n{'='*60}", flush=True)
    if has_critical:
        print("RESULT: critical issues found — resolve before merging", flush=True)
    else:
        print("RESULT: no critical issues", flush=True)
    print(f"{'='*60}\n", flush=True)

    return 1 if has_critical else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hermes workflow code review")
    parser.add_argument("workflows", nargs="*",
                        help="Workflow names to review (default: all)")
    parser.add_argument("--changed", action="store_true",
                        help="Auto-detect workflows changed in last git commit")
    args = parser.parse_args()

    if args.changed:
        targets = get_changed_workflows()
        if not targets:
            print("No workflow changes detected in last commit. Nothing to review.")
            sys.exit(0)
    elif args.workflows:
        targets = args.workflows
    else:
        targets = ALL_WORKFLOWS

    invalid = [t for t in targets if not (WF_DIR / t).exists()]
    if invalid:
        sys.exit(f"ERROR: Unknown workflow(s): {invalid}")

    sys.exit(asyncio.run(main(targets)))
