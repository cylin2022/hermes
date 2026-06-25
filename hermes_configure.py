#!/usr/bin/env python3
"""
hermes_configure.py — Platform auto-detection and config generation

Run once after installing Hermes on a new machine to:
  1. Detect CPU cores, RAM, GPUs
  2. Validate Docker + NVIDIA Container Toolkit permissions
  3. Generate an optimized config YAML for each workflow
  4. Print a summary of what was found and what needs attention

Usage:
    python3 hermes_configure.py [--workflow wgs_snp] [--output /path/to/config.yaml]
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERMES_DIR = Path(__file__).parent


# ── Hardware detection ─────────────────────────────────────────────────────────

def detect_cpu() -> int:
    try:
        return int(subprocess.check_output(["nproc"], text=True).strip())
    except Exception:
        return os.cpu_count() or 8


def detect_ram_gb() -> int:
    try:
        out = subprocess.check_output(["free", "-g"], text=True)
        for line in out.splitlines():
            if line.startswith("Mem:"):
                return int(line.split()[1])
    except Exception:
        pass
    return 64


def detect_gpus() -> list[dict]:
    """Return list of GPU info dicts; empty list if no GPU or nvidia-smi absent."""
    if not shutil.which("nvidia-smi"):
        return []
    try:
        out = subprocess.check_output(
            ["nvidia-smi",
             "--query-gpu=index,name,memory.total,driver_version",
             "--format=csv,noheader,nounits"],
            text=True
        ).strip()
        gpus = []
        for line in out.splitlines():
            if not line.strip():
                continue
            idx, name, mem, drv = [x.strip() for x in line.split(",", 3)]
            gpus.append({"index": idx, "name": name,
                          "memory_mib": int(mem), "driver": drv})
        return gpus
    except Exception:
        return []


# ── Permission checks ──────────────────────────────────────────────────────────

def check_docker() -> dict:
    result = {"installed": False, "user_in_group": False, "runnable": False, "error": ""}
    if not shutil.which("docker"):
        result["error"] = "docker not found in PATH"
        return result
    result["installed"] = True

    import grp
    try:
        docker_gid = grp.getgrnam("docker").gr_gid
        result["user_in_group"] = docker_gid in os.getgroups()
    except KeyError:
        result["error"] = "docker group does not exist"

    try:
        subprocess.check_output(
            ["docker", "info"], stderr=subprocess.DEVNULL, timeout=10
        )
        result["runnable"] = True
    except Exception as e:
        result["error"] = str(e)
    return result


def check_nvidia_container_toolkit(gpus: list) -> dict:
    result = {"available": False, "error": ""}
    if not gpus:
        result["error"] = "No GPU detected — NVIDIA toolkit not needed"
        return result
    try:
        subprocess.check_output(
            ["docker", "run", "--rm", "--gpus", "all",
             "nvidia/cuda:12.0-base-ubuntu22.04", "nvidia-smi"],
            stderr=subprocess.DEVNULL, timeout=60
        )
        result["available"] = True
    except Exception as e:
        result["error"] = f"GPU container test failed: {e}"
    return result


# ── Parameter recommendation ───────────────────────────────────────────────────

def recommend_params(cores: int, ram_gb: int, gpus: list) -> dict:
    """
    Recommend workflow parameters for this machine.

    threads: main CPU parallelism for BWA-MEM2, samtools, etc.
    dv_shards: DeepVariant --num_shards (CPU-intensive make_examples phase)
    mem_gb: memory per DeepVariant job

    Strategy:
    - Reserve ~10% of cores for OS / other processes
    - When GPU present: threads = floor(cores * 0.4) so 2 BWA-MEM2 + 1 DV
      can run simultaneously without oversubscription
    - When no GPU: threads = floor(cores * 0.5)
    - dv_shards = threads (DV makes_examples uses all shards in parallel)
    - mem_gb = min(256, ram_gb // 4)  (DV peaks at ~4 GB/shard × shards)
    """
    usable = int(cores * 0.9)

    if gpus:
        # GPU path: BWA-MEM2 runs on CPU, DV runs on GPU
        # We want: 2 × threads (BWA) + spare ≤ usable
        # threads = floor(usable / 2.5) → leaves room for 1 extra job
        threads = max(8, usable // 2)
        dv_shards = threads
    else:
        threads = max(8, usable // 2)
        dv_shards = threads

    # DV memory: 4 GB/shard during make_examples, but shards are sequential within a job
    # Peak = ~4 GB × num_shards / parallelism_factor ≈ 4 GB × shards × 0.6
    # Cap at 512 GB to be safe
    dv_mem = min(512, max(64, dv_shards * 4))
    dv_mem = min(dv_mem, ram_gb // 2)  # never more than half system RAM

    return {
        "threads": threads,
        "dv_shards": dv_shards,
        "mem_gb": dv_mem,
        "use_gpu": bool(gpus),
    }


# ── Config generation ─────────────────────────────────────────────────────────

def generate_wgs_snp_config(params: dict, output_path: Path | None = None) -> str:
    template_path = HERMES_DIR / "workflows" / "wgs_snp" / "config_template.yaml"
    template = template_path.read_text() if template_path.exists() else ""

    generated = f"""\
# Auto-generated by hermes_configure.py
# Machine: {detect_cpu()} cores, {detect_ram_gb()} GB RAM, {len(detect_gpus())} GPU(s)
# Edit paths below before running.

# ── Input ─────────────────────────────────────────────────────────────────────
samplesheet: "/path/to/samplesheet.csv"

# ── Reference ─────────────────────────────────────────────────────────────────
genome_fasta: "/path/to/genome.fa"
species_name: "my_species"

# ── Annotation ────────────────────────────────────────────────────────────────
gtf: "/path/to/annotation.gtf"
gff: ""

# ── Variant calling ───────────────────────────────────────────────────────────
ploidy: 2

# ── GPU / Docker ──────────────────────────────────────────────────────────────
use_gpu: {"true" if params["use_gpu"] else "false"}
dv_docker_image: "google/deepvariant:1.10.0-gpu"

# ── Resources (auto-tuned for this machine) ───────────────────────────────────
threads: {params["threads"]}          # CPU threads for BWA-MEM2, samtools
dv_shards: {params["dv_shards"]}      # DeepVariant --num_shards (independent of threads)
                                       # threads: 4 is used for Snakemake scheduling only
mem_gb: {params["mem_gb"]}            # memory per DeepVariant job (GB)

# ── Quality filters ───────────────────────────────────────────────────────────
min_gq: 20
min_depth: 5

# ── Output ────────────────────────────────────────────────────────────────────
outdir: "/path/to/output"
"""
    if output_path:
        output_path.write_text(generated)
        print(f"Config written to: {output_path}")
    return generated


# ── Report ────────────────────────────────────────────────────────────────────

def print_report(cores, ram_gb, gpus, docker, nct, params):
    OK  = "\033[92m✓\033[0m"
    WARN = "\033[93m⚠\033[0m"
    ERR = "\033[91m✗\033[0m"

    def status(ok, warn=False):
        return OK if ok else (WARN if warn else ERR)

    print("\n=== Hermes Platform Report ===\n")

    # Hardware
    print(f"CPU cores : {cores}")
    print(f"RAM       : {ram_gb} GB")
    if gpus:
        for g in gpus:
            print(f"GPU [{g['index']}]  : {g['name']}  {g['memory_mib']} MiB  driver {g['driver']}")
    else:
        print("GPU       : none detected")

    # Permissions
    print()
    print(f" {status(docker['installed'])} Docker installed")
    print(f" {status(docker['user_in_group'])} Current user in docker group")
    print(f" {status(docker['runnable'])} Docker daemon accessible")
    if not docker["user_in_group"]:
        print(f"    → Fix: sudo usermod -aG docker $USER && newgrp docker")
    if not docker["runnable"] and docker["error"]:
        print(f"    → Error: {docker['error']}")

    if gpus:
        if nct["error"] == "skipped":
            print(f" - NVIDIA Container Toolkit (skipped — rerun without --skip-gpu-test to verify)")
        else:
            print(f" {status(nct['available'])} NVIDIA Container Toolkit (GPU passthrough)")
            if not nct["available"]:
                print(f"    → Error: {nct['error']}")
                print(f"    → Fix: see INSTALL.md section 2.4")

    # Recommendations
    print()
    print("=== Recommended Parameters ===\n")
    print(f"  threads  : {params['threads']:<6}  (BWA-MEM2, samtools)")
    print(f"  dv_shards: {params['dv_shards']:<6}  (DeepVariant make_examples)")
    print(f"  mem_gb   : {params['mem_gb']:<6}  (DeepVariant memory cap)")
    print(f"  use_gpu  : {str(params['use_gpu']).lower()}")

    # Scheduling note
    if gpus:
        bwa_pairs = cores // params["threads"]
        print(f"\n  With these settings:")
        print(f"    Up to {bwa_pairs} BWA-MEM2 jobs × {params['threads']} threads")
        print(f"    + DeepVariant: threads:1 (scheduling), {params['dv_shards']} shards, 1 GPU")
        print(f"    Concurrency rule: BWA-MEM2 limited by CPU; DeepVariant limited by GPU")
        print(f"    → They run fully in parallel, no resource conflict")

    print()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Hermes platform auto-configuration")
    parser.add_argument("--workflow", default="wgs_snp",
                        choices=["wgs_snp"],
                        help="Workflow to configure")
    parser.add_argument("--output", "-o", default=None,
                        help="Write generated config to this path")
    parser.add_argument("--skip-gpu-test", action="store_true",
                        help="Skip the slow GPU container test")
    args = parser.parse_args()

    print("Detecting platform...", end=" ", flush=True)
    cores  = detect_cpu()
    ram_gb = detect_ram_gb()
    gpus   = detect_gpus()
    print("done")

    print("Checking Docker permissions...", end=" ", flush=True)
    docker = check_docker()
    print("done")

    nct = {"available": False, "error": "skipped"}
    if gpus and not args.skip_gpu_test and docker["runnable"]:
        print("Testing GPU container passthrough (may take ~30s)...", end=" ", flush=True)
        nct = check_nvidia_container_toolkit(gpus)
        print("done")

    params = recommend_params(cores, ram_gb, gpus)
    print_report(cores, ram_gb, gpus, docker, nct, params)

    output_path = Path(args.output) if args.output else None
    if args.workflow == "wgs_snp":
        config_text = generate_wgs_snp_config(params, output_path)
        if not output_path:
            print("=== Generated Config (wgs_snp) ===\n")
            print(config_text)


if __name__ == "__main__":
    main()
