#!/usr/bin/env python3
"""
Hermes MCP Server — Bioinformatics Pipeline Bridge
Exposes local workflow control to Claude Code via MCP stdio protocol.

Transport: stdio only (launched by Claude Code, not network-accessible).
"""

import asyncio
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

from fastmcp import FastMCP

HERMES_DIR = Path(__file__).parent
RUNS_DIR   = HERMES_DIR / "runs"
WF_DIR     = HERMES_DIR / "workflows"
RUNS_DIR.mkdir(exist_ok=True)

mcp = FastMCP("Hermes")

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


# ── Workflow execution ────────────────────────────────────────────────────────

@mcp.tool()
def list_workflows() -> list[dict]:
    """List available Snakemake workflow templates."""
    workflows = []
    for d in sorted(WF_DIR.iterdir()):
        if d.is_dir() and (d / "Snakefile").exists():
            readme = (d / "README.md").read_text() if (d / "README.md").exists() else ""
            workflows.append({
                "name":            d.name,
                "description":     readme.splitlines()[0].lstrip("# ") if readme else d.name,
                "config_template": str(d / "config_template.yaml")
                    if (d / "config_template.yaml").exists() else None,
            })
    return workflows


@mcp.tool()
def run_workflow(
    workflow_name: str,
    config: dict,
    run_id: str | None = None,
    cores: int = 128,
    dry_run: bool = False,
) -> dict:
    """
    Start a Snakemake workflow.

    Args:
        workflow_name: Name matching a directory under workflows/
        config: Key-value pairs passed as Snakemake --config
        run_id: Optional label (auto-generated if omitted)
        cores: CPU cores to allocate (default 128 of 160 available)
        dry_run: If True, preview steps without executing
    """
    snakefile = WF_DIR / workflow_name / "Snakefile"
    if not snakefile.exists():
        return {"error": f"Workflow '{workflow_name}' not found"}

    run_id  = run_id or f"{workflow_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    log_path  = run_dir / "snakemake.log"
    pid_path  = run_dir / "pid"
    meta_path = run_dir / "meta.json"

    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--cores",     str(cores),
        "--use-conda",
        "--rerun-incomplete",
        "--config", *[f"{k}={v}" for k, v in config.items()],
        *(["--dry-run"] if dry_run else []),
    ]

    if dry_run:
        result = subprocess.run(cmd, capture_output=True, text=True)
        return {"dry_run_output": result.stdout, "stderr": result.stderr}

    with open(log_path, "w") as log_file:
        proc = subprocess.Popen(
            cmd, stdout=log_file, stderr=subprocess.STDOUT,
            cwd=str(run_dir),
        )

    pid_path.write_text(str(proc.pid))
    (run_dir / "snakemake.pid").write_text(str(proc.pid))
    meta_path.write_text(json.dumps({
        "run_id":     run_id,
        "workflow":   workflow_name,
        "config":     config,
        "cores":      cores,
        "started_at": datetime.now().isoformat(),
        "pid":        proc.pid,
    }, indent=2))

    return {"run_id": run_id, "pid": proc.pid, "log": str(log_path), "status": "started"}


# ── Monitoring ────────────────────────────────────────────────────────────────

@mcp.tool()
def list_runs(limit: int = 20) -> list[dict]:
    """List recent workflow runs with status."""
    runs = []
    for run_dir in sorted(RUNS_DIR.iterdir(), reverse=True)[:limit]:
        if not run_dir.is_dir():
            continue
        meta_path = run_dir / "meta.json"
        if not meta_path.exists():
            continue
        meta = json.loads(meta_path.read_text())
        runs.append({**meta, "status": _run_status(run_dir / "pid", run_dir)})
    return runs


@mcp.tool()
def get_status(run_id: str) -> dict:
    """Get status and recent log tail for a workflow run."""
    run_dir = RUNS_DIR / run_id
    if not (run_dir / "meta.json").exists():
        return {"error": f"Run '{run_id}' not found"}

    meta     = json.loads((run_dir / "meta.json").read_text())
    status   = _run_status(run_dir / "pid", run_dir)
    log_path = run_dir / "snakemake.log"
    log_tail = "\n".join(log_path.read_text().splitlines()[-30:]) if log_path.exists() else ""

    return {**meta, "status": status, "log_tail": log_tail}


@mcp.tool()
def get_log(run_id: str, lines: int = 50) -> str:
    """Get the last N lines of a run's log."""
    log_path = RUNS_DIR / run_id / "snakemake.log"
    if not log_path.exists():
        return f"No log found for run '{run_id}'"
    return "\n".join(log_path.read_text().splitlines()[-lines:])


@mcp.tool()
def stop_run(run_id: str) -> dict:
    """Stop a running workflow (SIGTERM to Snakemake process)."""
    pid_path = RUNS_DIR / run_id / "pid"
    if not pid_path.exists():
        return {"error": f"Run '{run_id}' not found or already finished"}
    pid = int(pid_path.read_text().strip())
    try:
        subprocess.run(["kill", "-TERM", str(pid)], check=False)
        return {"run_id": run_id, "pid": pid, "status": "stop signal sent"}
    except Exception as e:
        return {"error": str(e)}


# ── Filesystem helpers ────────────────────────────────────────────────────────

@mcp.tool()
def read_file(path: str, tail_lines: int = 100) -> str:
    """Read the last N lines of a file (logs, result tables, etc.)."""
    p = Path(path)
    if not p.exists():
        return f"File not found: {path}"
    content = p.read_text(errors="replace").splitlines()
    return "\n".join(content[-tail_lines:])


@mcp.tool()
def list_files(directory: str, pattern: str = "*") -> list[str]:
    """List files matching a glob pattern in a directory."""
    d = Path(directory)
    if not d.exists():
        return [f"Directory not found: {directory}"]
    return sorted(str(p) for p in d.glob(pattern))


# ── Workflow critique (LLM-based code review) ────────────────────────────────

def _bundle_workflow_files(wf_dir: Path) -> str:
    parts = []
    candidates = [
        (wf_dir / "Snakefile",             None),
        (wf_dir / "config_template.yaml",  None),
        (wf_dir / "envs",                  "*.yaml"),
        (wf_dir / "scripts",               "*.R"),
        (wf_dir / "scripts",               "*.py"),
    ]
    for target, glob_pat in candidates:
        if glob_pat:
            if target.exists():
                for f in sorted(target.glob(glob_pat)):
                    parts.append(f"=== {f.relative_to(wf_dir)} ===\n{f.read_text(errors='replace')}")
        elif target.exists():
            parts.append(f"=== {target.name} ===\n{target.read_text(errors='replace')}")
    return "\n\n".join(parts)


def _critique_prompt(name: str, files_content: str, checklist: str) -> str:
    return f"""You are a senior bioinformatics engineer doing an independent code review of the
"{name}" Snakemake workflow. Review with fresh eyes — no knowledge of design decisions.
Find bugs that cause crashes or wrong results in production.
Cost of a missed bug: days of lost compute on a 160-core / 2.2 TiB server.

## Workflow files

{files_content}

## PIPELINE_CHECKLIST.md

{checklist}

## What to check

**Docker (CRITICAL):** --user $(id -u):$(id -g) on every docker run (EXCEPTION: BRAKER3 root).
mkdir -p output dir BEFORE docker run. busybox chown for root-owned outputs.

**Shell blocks (CRITICAL):** set -euo pipefail first line. test -s {{output}} after every command.
sed back-references: use \\\\1 not \\1 (Python processes Snakemake shell strings first).

**samtools order (CRITICAL):** fixmate requires name-sorted input (sort -n before fixmate).

**Index files:** .fai/.bai/.csi/.tbi declared as inputs where tools require them.
samtools >= 1.12 writes .csi (not .bai) for non-standard references.

**Conda envs:** every library()/import must have a matching package in the conda yaml.

**Statistics:** appropriate normalization, multiple-testing correction, NA/empty-output guards.

Use the submit_review tool. Classify CRITICAL (crash or wrong result) vs MINOR (style/sub-optimal).
Include file:line for every issue."""


async def _critique_one(client, name: str, checklist: str) -> dict | None:
    wf_dir = WF_DIR / name
    if not wf_dir.exists():
        return None
    files_content = _bundle_workflow_files(wf_dir)

    import anthropic as _anthropic
    try:
        response = await client.messages.create(
            model       = "claude-sonnet-4-6",
            max_tokens  = 4096,
            tools       = [{
                "name":         "submit_review",
                "description":  "Submit structured code review findings",
                "input_schema": ISSUE_SCHEMA,
            }],
            tool_choice = {"type": "tool", "name": "submit_review"},
            messages    = [{"role": "user",
                            "content": _critique_prompt(name, files_content, checklist)}],
        )
    except _anthropic.APIError as exc:
        return {"workflow": name, "critical_issues": [], "minor_issues": [],
                "verdict": "minor_fixes_needed", "summary": f"API error: {exc}"}

    for block in response.content:
        if block.type == "tool_use" and block.name == "submit_review":
            result = dict(block.input)
            result["workflow"] = name
            return result
    return None


async def _synthesize(client, results: list[dict]) -> str:
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
            f"CRITICAL ({len(r.get('critical_issues',[]))}):\n{crit}\n"
            f"MINOR ({len(r.get('minor_issues',[]))}):\n{minor}\n"
            f"Summary: {r.get('summary','')}"
        )

    response = await client.messages.create(
        model      = "claude-sonnet-4-6",
        max_tokens = 2048,
        messages   = [{"role": "user", "content": (
            "Synthesize these code-review results for Hermes bioinformatics workflows.\n\n"
            + "\n\n".join(details)
            + "\n\n"
            "## 1. Overall verdict (1-2 sentences)\n"
            "## 2. Workflows ready to run\n"
            "## 3. Critical fixes required (group by workflow; exact code fixes)\n"
            "## 4. Systemic patterns (issues in 2+ workflows)\n"
            "## 5. Next steps\n"
            "Be specific. Exact code where possible. No filler."
        )}],
    )
    return response.content[0].text


@mcp.tool()
async def critique_workflow(workflows: list[str] | None = None) -> dict:
    """
    Run independent LLM code review of Hermes workflow(s) before production.

    Equivalent to: Workflow({ scriptPath: critique.js, args: workflows })
    Requires ANTHROPIC_API_KEY in the MCP server environment.

    Args:
        workflows: Workflow names to review, e.g. ["pool_seq", "rnaseq"].
                   Omit to review all 11 workflows (~2-4 min).
    Returns:
        { findings, synthesis, stats: {n_reviewed, n_critical, n_minor} }
    """
    import anthropic as _anthropic

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return {"error": "ANTHROPIC_API_KEY not set in MCP server environment"}

    targets = workflows if workflows else ALL_WORKFLOWS
    invalid = [w for w in targets if not (WF_DIR / w).exists()]
    if invalid:
        return {"error": f"Unknown workflow(s): {invalid}. Available: {ALL_WORKFLOWS}"}

    checklist = (HERMES_DIR / "PIPELINE_CHECKLIST.md").read_text() \
        if (HERMES_DIR / "PIPELINE_CHECKLIST.md").exists() else ""

    client  = _anthropic.AsyncAnthropic(api_key=api_key)
    results = await asyncio.gather(
        *[_critique_one(client, name, checklist) for name in targets],
        return_exceptions=True,
    )

    valid      = [r for r in results if isinstance(r, dict)]
    n_critical = sum(len(r.get("critical_issues", [])) for r in valid)
    n_minor    = sum(len(r.get("minor_issues",    [])) for r in valid)

    synthesis = await _synthesize(client, valid) if valid else ""

    return {
        "findings":  valid,
        "synthesis": synthesis,
        "stats":     {"n_reviewed": len(valid), "n_critical": n_critical, "n_minor": n_minor},
    }


# ── Internal helpers ─────────────────────────────────────────────────────────

def _run_status(pid_path: Path, run_dir: Path) -> str:
    if not pid_path.exists():
        return "unknown"
    pid = int(pid_path.read_text().strip())
    if Path(f"/proc/{pid}").exists():
        return "running"
    log_path = run_dir / "snakemake.log"
    if log_path.exists():
        tail = log_path.read_text().splitlines()[-5:]
        if any("(100%) done" in l or "steps (100%) done" in l for l in tail):
            return "completed"
        if any("Error" in l or "error" in l or "ERROR" in l for l in tail):
            return "failed"
    return "finished"


# ── Entry point — stdio only ─────────────────────────────────────────────────

if __name__ == "__main__":
    mcp.run(transport="stdio")
