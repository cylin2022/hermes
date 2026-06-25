#!/usr/bin/env python3
"""
Hermes MCP Server — Bioinformatics Pipeline Bridge
Exposes local workflow control to Claude Code via MCP stdio protocol.

Transport: stdio only (launched by Claude Code, not network-accessible).
"""

import json
import subprocess
from datetime import datetime
from pathlib import Path

from fastmcp import FastMCP

HERMES_DIR = Path(__file__).parent
RUNS_DIR   = HERMES_DIR / "runs"
WF_DIR     = HERMES_DIR / "workflows"
RUNS_DIR.mkdir(exist_ok=True)

mcp = FastMCP("Hermes")


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
