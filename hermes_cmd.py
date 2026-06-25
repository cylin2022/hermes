#!/usr/bin/env python3
"""
小賀 指令監聽 daemon
訂閱 ntfy 指令頻道，執行預定義指令，結果推回通知頻道。

支援指令（手機傳送純文字）：
  status / 狀態     — 目前流程進度 + 資源狀況
  log / 日誌        — Snakemake log 最後 30 行
  gpu               — GPU 使用狀況
  cpu               — CPU / 記憶體概況
  jobs              — 目前跑中的 job 數
  stop / 停止       — 送 SIGTERM 給 Snakemake（需再確認）
  stop confirm      — 確認停止
  help / 說明       — 顯示指令列表
"""
import json
import os
import re
import signal
import subprocess
import time
import urllib.request
from pathlib import Path

# ── config ────────────────────────────────────────────────────────────────────
CFG_FILE    = Path.home() / ".hermes_notify.json"
cfg         = json.loads(CFG_FILE.read_text())
NTFY_BASE   = "https://ntfy.sh"
CMD_TOPIC   = cfg.get("ntfy_cmd_topic",    "ntfy_cylin-hermes2026-cmd")
CMD_PREFIX  = cfg.get("cmd_password",      "yama2026")
NOTIFY_TOPIC = cfg.get("ntfy_topic",       "ntfy_cylin-hermes2026")
LOG_FILE    = Path(cfg.get("wgs_log",
    "/home/cylin/hermes/runs/wgs_snp_tilapia_20260613/snakemake.log"))
PID_FILE    = Path(cfg.get("wgs_pid_file",
    "/home/cylin/hermes/runs/wgs_snp_tilapia_20260613/snakemake.pid"))

pending_stop = False   # waiting for "stop confirm"


# ── notify helpers ────────────────────────────────────────────────────────────
def push(title: str, msg: str, tags: str = "robot", priority: str = "default") -> None:
    try:
        subprocess.run(
            ["curl", "-s",
             "-H", f"Title: {title}",
             "-H", f"Tags: {tags}",
             "-H", f"Priority: {priority}",
             "-d", msg,
             f"{NTFY_BASE}/{NOTIFY_TOPIC}"],
            check=False, capture_output=True
        )
    except Exception as e:
        print(f"[push error] {e}", flush=True)


# ── snakemake helpers ─────────────────────────────────────────────────────────
def find_snakemake_pid() -> int | None:
    if PID_FILE.exists():
        try:
            pid = int(PID_FILE.read_text().strip())
            os.kill(pid, 0)   # verify alive; raises if stale
            return pid
        except (ValueError, ProcessLookupError, PermissionError):
            pass   # stale PID — fall through to pgrep
    # fallback: find by process name
    r = subprocess.run(
        ["pgrep", "-f", "snakemake.*wgs_snp"],
        capture_output=True, text=True
    )
    pids = r.stdout.strip().splitlines()
    return int(pids[0]) if pids else None


def snakemake_running() -> bool:
    pid = find_snakemake_pid()
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False


def log_tail(n: int = 30) -> str:
    if not LOG_FILE.exists():
        return "（log 尚不存在）"
    lines = LOG_FILE.read_text().splitlines()
    return "\n".join(lines[-n:])


def count_finished(rule: str) -> int:
    if not LOG_FILE.exists():
        return 0
    text = LOG_FILE.read_text()
    # Count "Finished job N" after each "rule <name>:" block
    return len(re.findall(rf"rule {rule}:", text))


def get_status() -> str:
    running = snakemake_running()
    status  = "🟢 執行中" if running else "🔴 已停止"

    # Count key rules from log
    lines   = []
    for rule, total, label in [
        ("fastp",           60, "fastp QC"),
        ("bwa_mem2",        60, "BWA-MEM2"),
        ("sort_markdup",    60, "sort/markdup"),
        ("deepvariant",     60, "DeepVariant"),
        ("glnexus",          1, "GLnexus"),
        ("snpeff_annotate",  1, "SNPeff"),
    ]:
        n = count_finished(rule)
        bar = "✅" if n >= total else f"{n}/{total}"
        lines.append(f"  {label}: {bar}")

    # GPU
    gpu = subprocess.run(
        ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.free",
         "--format=csv,noheader,nounits"],
        capture_output=True, text=True
    ).stdout.strip()

    return (
        f"Snakemake: {status}\n\n"
        + "\n".join(lines)
        + f"\n\nGPU: {gpu} (util%, used MiB, free MiB)"
    )


def get_resources() -> str:
    gpu = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,utilization.gpu,memory.used,memory.free",
         "--format=csv,noheader"],
        capture_output=True, text=True
    ).stdout.strip()

    mem = subprocess.run(
        ["free", "-h"], capture_output=True, text=True
    ).stdout.strip()

    load = subprocess.run(
        ["uptime"], capture_output=True, text=True
    ).stdout.strip()

    return f"GPU:\n{gpu}\n\nRAM:\n{mem}\n\n{load}"


# ── command handler ───────────────────────────────────────────────────────────
def handle(msg: str) -> None:
    global pending_stop
    msg = msg.strip()

    # password prefix check
    if not msg.lower().startswith(CMD_PREFIX.lower()):
        print(f"[cmd] rejected (no prefix): {repr(msg)}", flush=True)
        return
    msg = msg[len(CMD_PREFIX):].strip()

    cmd = msg.lower()
    print(f"[cmd] {repr(cmd)}", flush=True)

    if cmd in ("help", "說明", "?", "？"):
        push("支援指令", (
            "status / 狀態  — 流程進度\n"
            "log / 日誌     — log 最後 30 行\n"
            "gpu            — GPU 狀況\n"
            "cpu            — CPU/RAM\n"
            "jobs           — 跑中的 job\n"
            "stop / 停止    — 停止流程（需再輸入 stop confirm）\n"
            "help / 說明    — 本說明"
        ), tags="memo")

    elif cmd in ("status", "狀態", "進度"):
        push("流程狀態", get_status(), tags="bar_chart")

    elif cmd in ("log", "日誌", "log30"):
        n = 30 if "30" not in cmd else 30
        push("Log 最後 30 行", log_tail(30), tags="scroll")

    elif cmd == "log50":
        push("Log 最後 50 行", log_tail(50), tags="scroll")

    elif cmd in ("gpu",):
        r = subprocess.run(
            ["nvidia-smi"], capture_output=True, text=True
        ).stdout.strip()
        push("GPU 狀況", r, tags="computer")

    elif cmd in ("cpu", "資源", "resource"):
        push("系統資源", get_resources(), tags="computer")

    elif cmd in ("jobs",):
        r = subprocess.run(
            ["ps", "aux"], capture_output=True, text=True
        ).stdout
        snk_lines = [l for l in r.splitlines() if "snakemake" in l or "deepvariant" in l.lower()]
        push("目前 Jobs", "\n".join(snk_lines[:20]) or "無 snakemake 程序", tags="mag")

    elif cmd in ("stop", "停止"):
        pending_stop = True
        push("⚠️ 確認停止？",
             "再傳送「stop confirm」確認停止 Snakemake。\n傳送其他任何指令取消。",
             tags="warning", priority="high")

    elif cmd == "stop confirm":
        if pending_stop:
            pending_stop = False
            pid = find_snakemake_pid()
            if pid:
                os.kill(pid, signal.SIGTERM)
                push("Snakemake 已停止", f"已送出 SIGTERM (PID {pid})", tags="octagonal_sign")
            else:
                push("找不到 Snakemake", "程序可能已結束", tags="question")
        else:
            push("請先輸入 stop", "需先輸入 stop 再輸入 stop confirm", tags="x")

    else:
        pending_stop = False   # any other command cancels stop
        push("未知指令", f"未知：{msg}\n輸入 help 查看支援指令", tags="x")


# ── main loop: subscribe to command topic ─────────────────────────────────────
def main():
    push("小賀 指令監聽啟動",
         f"指令頻道：{CMD_TOPIC}\n輸入 help 查看可用指令",
         tags="satellite")

    url = f"{NTFY_BASE}/{CMD_TOPIC}/json"
    print(f"[cmd daemon] subscribing to {url}", flush=True)

    while True:
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/x-ndjson"})
            with urllib.request.urlopen(req, timeout=300) as resp:
                for raw in resp:
                    try:
                        event = json.loads(raw.decode("utf-8"))
                    except json.JSONDecodeError:
                        continue
                    if event.get("event") == "message":
                        handle(event.get("message", ""))
        except Exception as e:
            print(f"[cmd daemon] connection error: {e} — retrying in 15s", flush=True)
            time.sleep(15)


if __name__ == "__main__":
    main()
