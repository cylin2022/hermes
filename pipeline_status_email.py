#!/usr/bin/env python3
"""
Collect status of running pipelines and send a summary email.
Called every 8 hours by Claude Code cron.
"""
import subprocess
import re
from datetime import datetime
from pathlib import Path

# ── Helpers ──────────────────────────────────────────────────────────────────

def run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return ""

def tail(path, n=10):
    return run(f"tail -{n} '{path}' 2>/dev/null")

# ── WGS SNP status ────────────────────────────────────────────────────────────

def wgs_status():
    gvcf_dir = Path("/home/cylin/WGS_AntiSalt_tilapia/wgs_snp_out/gvcf")
    done = sorted(gvcf_dir.glob("*.g.vcf.gz")) if gvcf_dir.exists() else []
    total = 60

    # Current DeepVariant sample from ps
    ps = run("ps aux | grep 'run_deepvariant' | grep -v grep | grep 'output_gvcf'")
    current_sample = ""
    m = re.search(r"output_gvcf=\S+/(\S+?)\.g\.vcf\.gz", ps)
    if m:
        current_sample = m.group(1)

    # Estimate time remaining based on last 4 completed samples
    recent = sorted(done, key=lambda p: p.stat().st_mtime)[-4:]
    avg_sec = None
    if len(recent) >= 2:
        intervals = []
        for i in range(1, len(recent)):
            intervals.append(recent[i].stat().st_mtime - recent[i-1].stat().st_mtime)
        avg_sec = sum(intervals) / len(intervals)

    remaining = total - len(done)
    eta_str = "—"
    if avg_sec and remaining > 0:
        eta_hours = (remaining * avg_sec) / 3600
        eta_str = f"~{eta_hours:.1f} 小時"

    lines = [
        f"【WGS SNP — wgs_snp_tilapia_20260613】",
        f"  完成：{len(done)} / {total} 個 gVCF",
        f"  目前處理：{current_sample or '（無 DeepVariant 程序）'}",
        f"  預估剩餘：{eta_str}",
    ]
    if done:
        last = done[-1] if not current_sample else sorted(done, key=lambda p: p.stat().st_mtime)[-1]
        lines.append(f"  最新完成：{last.stem.replace('.g.vcf', '')}")

    running = bool(current_sample)
    return "\n".join(lines), running


# ── Genome Annotation (BRAKER3) status ────────────────────────────────────────

def braker_status():
    log = Path("/home/cylin/TWN_Hamaguri/assembly_output/annotation/braker4/braker.log")
    gtf = Path("/home/cylin/TWN_Hamaguri/assembly_output/annotation/braker4/braker.gtf")

    if gtf.exists() and gtf.stat().st_size > 0:
        return "【基因組注釋 — TWN_Hamaguri】\n  ✅ BRAKER4 已完成，braker.gtf 產出", True  # completed = still "ok", not stopped

    # Check if braker3 container is running
    ps = run("ps aux | grep 'braker.pl' | grep -v grep")
    running = bool(ps)

    stage = "—"
    if log.exists():
        recent = tail(str(log), 20)
        # Find last meaningful step
        for line in reversed(recent.splitlines()):
            line = line.strip()
            if line.startswith("#") and ":" in line:
                stage = line.split("2026:")[-1].strip() if "2026:" in line else line
                break

    # Start time from log
    start_str = ""
    if log.exists():
        first = run(f"head -5 '{log}'")
        m = re.search(r"Jun \d+ (\d+:\d+:\d+) 2026", first)
        if m:
            start_str = m.group(1)

    elapsed = ""
    if start_str:
        try:
            t0 = datetime.strptime(f"2026-06-18 {start_str}", "%Y-%m-%d %H:%M:%S")
            diff = (datetime.now() - t0).total_seconds() / 3600
            elapsed = f"{diff:.1f} 小時"
        except Exception:
            pass

    lines = [
        f"【基因組注釋 — TWN_Hamaguri】",
        f"  狀態：{'執行中' if running else '⚠️ 無 BRAKER3 程序'}",
        f"  已執行：{elapsed or '—'}",
        f"  目前步驟：{stage}",
        f"  braker.gtf：{'尚未產出' if not gtf.exists() else '存在 (' + str(gtf.stat().st_size // 1024) + ' KB)'}",
    ]
    return "\n".join(lines), running


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    now = datetime.now().strftime("%Y-%m-%d %H:%M CST")

    wgs_report, wgs_running = wgs_status()
    braker_report, braker_running = braker_status()

    # braker_running=True means "completed OK" (not a running process).
    # wgs_running=True means DeepVariant is actively running.
    # both_done = braker completed AND wgs has no more gVCFs to process.
    gvcf_done = len(list(Path("/home/cylin/WGS_AntiSalt_tilapia/wgs_snp_out/gvcf").glob("*.g.vcf.gz"))) if Path("/home/cylin/WGS_AntiSalt_tilapia/wgs_snp_out/gvcf").exists() else 0
    wgs_complete = gvcf_done >= 60 and not wgs_running
    braker_complete = braker_running  # braker_running=True means gtf exists = completed
    both_done = wgs_complete and braker_complete
    # Only warn if a pipeline appears unexpectedly stopped (not simply completed)
    wgs_unexpectedly_stopped = not wgs_running and gvcf_done < 60
    subject = f"分析進度報告 {now}" if not both_done else f"✅ 所有分析已完成 {now}"

    body = f"""小賀 定期進度報告
時間：{now}
{'=' * 48}

{wgs_report}

{braker_report}

{'=' * 48}
{'⚠️  注意：WGS SNP 無 DeepVariant 程序，尚有 gVCF 未完成，請確認。' if wgs_unexpectedly_stopped else ''}
{'✅ 兩個流程均已完成。' if both_done else ''}
"""

    # Send via notify.py
    import sys
    notify = Path("/home/cylin/hermes/notify.py")
    result = subprocess.run(
        [sys.executable, str(notify), subject, body],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        print(f"Email sent: {subject}")
    else:
        print(f"Email failed: {result.stderr}", file=sys.stderr)
        sys.exit(1)

    # If both done, signal completion
    if both_done:
        print("Both pipelines complete. Cron should be stopped.")


if __name__ == "__main__":
    main()
