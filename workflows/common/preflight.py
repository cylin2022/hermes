"""
Hermes preflight resource checker.

Call preflight(config, needs_gpu=False, needs_docker=False) inside
the Snakemake onstart: block of every workflow.

Checks performed:
  CPU    — requested threads vs available cores
  RAM    — requested mem_gb vs available RAM
  Disk   — output directory has sufficient free space
  GPU    — nvidia-smi reachable (when needs_gpu=True)
  Docker — daemon running (when needs_docker=True)

Policy:
  WARN  (yellow)  — resource is tight but pipeline can still attempt to run
  ERROR (red)     — resource is critically insufficient; pipeline is aborted
"""

import os
import shutil
import subprocess
import sys


# ── ANSI colours (suppressed in non-TTY environments) ─────────────────────────
def _c(code, text):
    return f"\033[{code}m{text}\033[0m" if sys.stderr.isatty() else text

RED    = lambda t: _c("1;31", t)
YELLOW = lambda t: _c("1;33", t)
GREEN  = lambda t: _c("1;32", t)
BOLD   = lambda t: _c("1;37", t)


def _section(title):
    print(f"\n{BOLD('━' * 60)}", file=sys.stderr)
    print(f"{BOLD(f'  Hermes preflight: {title}')}", file=sys.stderr)
    print(BOLD('━' * 60), file=sys.stderr)


def _ok(msg):    print(f"  {GREEN('✓')}  {msg}", file=sys.stderr)
def _warn(msg):  print(f"  {YELLOW('⚠')}  {YELLOW(msg)}", file=sys.stderr)
def _err(msg):   print(f"  {RED('✗')}  {RED(msg)}", file=sys.stderr)
def _info(msg):  print(f"     {msg}", file=sys.stderr)


# ── Individual checks ──────────────────────────────────────────────────────────

def check_cpu(requested_threads):
    available = os.cpu_count() or 1
    pct = requested_threads / available * 100
    _info(f"CPU cores available : {available}")
    _info(f"Threads requested   : {requested_threads}  ({pct:.0f}% of available)")
    if requested_threads > available:
        _warn(f"threads ({requested_threads}) exceeds available cores ({available}). "
              "Snakemake will over-subscribe — jobs may be slower.")
    elif pct >= 90:
        _warn(f"Using {pct:.0f}% of CPU. Other processes on this machine will be starved.")
    else:
        _ok(f"CPU OK ({requested_threads}/{available} cores)")
    return True   # CPU over-subscription is a warn, never an abort


def check_ram(requested_gb):
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    avail_kb = int(line.split()[1])
                    break
        avail_gb = avail_kb / 1024 / 1024
    except Exception:
        _warn("Cannot read /proc/meminfo — RAM check skipped")
        return True

    pct = requested_gb / avail_gb * 100
    _info(f"RAM available       : {avail_gb:.0f} GB")
    _info(f"mem_gb requested    : {requested_gb} GB  ({pct:.0f}% of available)")

    if requested_gb > avail_gb * 0.95:
        _err(f"mem_gb ({requested_gb} GB) exceeds 95% of available RAM ({avail_gb:.0f} GB). "
             "Pipeline will likely be killed by OOM killer.")
        return False
    elif requested_gb > avail_gb * 0.80:
        _warn(f"mem_gb ({requested_gb} GB) is {pct:.0f}% of available RAM. "
              "May trigger OOM if other processes are running.")
    else:
        _ok(f"RAM OK ({requested_gb} GB requested / {avail_gb:.0f} GB available)")
    return True


def check_disk(outdir, min_free_gb=200):
    try:
        os.makedirs(outdir, exist_ok=True)
        usage = shutil.disk_usage(outdir)
        free_gb = usage.free / 1024 ** 3
        total_gb = usage.total / 1024 ** 3
    except Exception as e:
        _warn(f"Cannot check disk for {outdir}: {e}")
        return True

    _info(f"Output dir          : {outdir}")
    _info(f"Disk free           : {free_gb:.0f} GB / {total_gb:.0f} GB")

    if free_gb < min_free_gb * 0.5:
        _err(f"Only {free_gb:.0f} GB free in {outdir}. "
             f"Minimum recommended: {min_free_gb} GB. Pipeline may fail mid-run.")
        return False
    elif free_gb < min_free_gb:
        _warn(f"Only {free_gb:.0f} GB free — less than recommended {min_free_gb} GB. "
              "Monitor disk usage during the run.")
    else:
        _ok(f"Disk OK ({free_gb:.0f} GB free)")
    return True


def check_gpu():
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total,memory.free",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr)
        for i, line in enumerate(result.stdout.strip().splitlines()):
            name, total, free = [x.strip() for x in line.split(",")]
            _info(f"GPU {i}               : {name}  total={total}  free={free}")
        _ok(f"GPU OK ({len(result.stdout.strip().splitlines())} device(s) found)")
        return True
    except FileNotFoundError:
        _err("nvidia-smi not found. GPU workflows require NVIDIA drivers.")
        return False
    except Exception as e:
        _err(f"GPU check failed: {e}")
        return False


def check_docker():
    try:
        result = subprocess.run(
            ["docker", "info", "--format", "{{.ServerVersion}}"],
            capture_output=True, text=True, timeout=15
        )
        if result.returncode != 0:
            raise RuntimeError("Docker daemon not responding")
        _ok(f"Docker OK (server v{result.stdout.strip()})")
        return True
    except FileNotFoundError:
        _err("docker not found in PATH.")
        return False
    except Exception as e:
        _err(f"Docker check failed: {e}")
        return False


# ── Main entry point ───────────────────────────────────────────────────────────

def preflight(config, workflow_name="workflow", needs_gpu=False, needs_docker=False):
    """
    Run all preflight checks. Call from Snakemake's onstart: block.
    Raises SystemExit(1) if any critical check fails.

    Parameters
    ----------
    config       : Snakemake config dict
    workflow_name: display name shown in the header
    needs_gpu    : set True for workflows that use DeepVariant GPU
    needs_docker : set True for workflows that invoke docker run
    """
    _section(workflow_name)

    threads     = int(config.get("threads", 64))
    mem_gb      = int(config.get("mem_gb", 256))
    outdir      = config.get("outdir", ".")
    min_disk_gb = int(config.get("min_disk_gb", 200))

    failures = []

    if not check_cpu(threads):
        failures.append("CPU")

    if not check_ram(mem_gb):
        failures.append("RAM")

    if not check_disk(outdir, min_free_gb=min_disk_gb):
        failures.append("DISK")

    if needs_gpu:
        if not check_gpu():
            failures.append("GPU")

    if needs_docker:
        if not check_docker():
            failures.append("DOCKER")

    print(file=sys.stderr)
    if failures:
        print(RED(f"  PREFLIGHT FAILED: {', '.join(failures)}"), file=sys.stderr)
        print(RED(  "  Fix the issues above before re-running."), file=sys.stderr)
        print(BOLD('━' * 60), file=sys.stderr)
        raise SystemExit(1)
    else:
        print(GREEN(f"  All preflight checks passed. Starting {workflow_name} …"),
              file=sys.stderr)
        print(BOLD('━' * 60), file=sys.stderr)
