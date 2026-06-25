#!/usr/bin/env python3
"""
Watch a wgs_snp Snakemake log and send notifications at key milestones.
Supports LINE Notify and/or Gmail. Config: ~/.hermes_notify.json
Usage: python3 monitor_wgs.py <log_file> [n_samples]
"""
import json
import re
import smtplib
import ssl
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

CONFIG    = Path.home() / ".hermes_notify.json"
LOG_FILE  = Path(sys.argv[1])
N_SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 60


def load_cfg():
    return json.loads(CONFIG.read_text()) if CONFIG.exists() else {}


def send_ntfy(title: str, msg: str, tags: str = "bell") -> None:
    cfg = load_cfg()
    topic = cfg.get("ntfy_topic")
    if not topic:
        return
    subprocess.run(
        ["curl", "-s",
         "-H", f"Title: {title}",
         "-H", f"Tags: {tags}",
         "-H", "Priority: high",
         "-d", msg,
         f"https://ntfy.sh/{topic}"],
        check=False
    )


def send_email(subject: str, body: str) -> None:
    cfg = load_cfg()
    if not cfg.get("gmail_user"):
        return
    msg_obj = __import__("email.mime.text", fromlist=["MIMEText"]).MIMEText(body, "plain", "utf-8")
    msg_obj["Subject"] = f"[小賀] {subject}"
    msg_obj["From"]    = cfg["gmail_user"]
    msg_obj["To"]      = cfg["notify_to"]
    ctx = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=ctx) as smtp:
        smtp.login(cfg["gmail_user"], cfg["gmail_app_password"])
        smtp.sendmail(cfg["gmail_user"], cfg["notify_to"], msg_obj.as_string())


def notify(subject: str, body: str, tags: str = "bell") -> None:
    try:
        send_ntfy(subject, body, tags)
    except Exception as e:
        print(f"ntfy failed: {e}", flush=True)
    try:
        send_email(subject, body)
    except Exception as e:
        print(f"Email notify failed: {e}", flush=True)


# ── milestones ────────────────────────────────────────────────────────────────
MILESTONES = {
    "fastp":           (N_SAMPLES, "QC/Trimming 完成",          "開始 BWA-MEM2 比對"),
    "bwa_mem2":        (N_SAMPLES, "BWA-MEM2 比對完成",          "開始 DeepVariant"),
    "sort_markdup":    (N_SAMPLES, "排序去重完成",               "DeepVariant 進行中"),
    "deepvariant":     (N_SAMPLES, "DeepVariant 完成",           "開始 GLnexus 聯合定型"),
    "glnexus":         (1,         "GLnexus 聯合定型完成",       "開始 SNPeff 注釋"),
    "snpeff_annotate": (1,         "SNPeff 注釋完成",            "流程即將結束"),
}

ERROR_RE  = re.compile(
    r"(Error in rule|CalledProcessError|"
    r"Exiting because a job execution failed|Killed|Out of memory|oom-kill)",
    re.I
)
# MissingOutputException excluded: Snakemake emits it during graceful shutdown / --rerun-incomplete
# FAILED excluded: too broad, matches normal bwa-mem2 orientation-skip messages

notified     = set()
counts       = defaultdict(int)
current_rule = ""
last_error   = 0.0   # timestamp of last error notify


def tail_forever(path: Path):
    while not path.exists():
        time.sleep(5)
    with open(path) as f:
        f.seek(0, 2)
        while True:
            line = f.readline()
            if line:
                yield line.rstrip("\n")
            else:
                time.sleep(10)


print(f"[monitor] Watching {LOG_FILE}  n_samples={N_SAMPLES}", flush=True)

for line in tail_forever(LOG_FILE):

    # Track which rule is currently being described
    m = re.search(r"^rule (\w+):", line)
    if m:
        current_rule = m.group(1)

    # Snakemake 9 prints "Finished job N." on its own line
    if "Finished job" in line and current_rule:
        counts[current_rule] += 1
        rule = current_rule
        if rule in MILESTONES and rule not in notified:
            expected, label, hint = MILESTONES[rule]
            if counts[rule] >= expected:
                notified.add(rule)
                body = (f"規則：{rule}  完成 {counts[rule]}/{expected}\n"
                        f"下一步：{hint}")
                notify(label, body, tags="white_check_mark")
                print(f"[notify] {label}", flush=True)

    # Workflow 100% done
    if re.search(r"\d+ of \d+ steps \(100%\) done", line):
        stats_f = Path("/home/cylin/WGS_AntiSalt_tilapia/wgs_snp_out/variant_stats.txt")
        stats   = stats_f.read_text().strip() if stats_f.exists() else "結果目錄：wgs_snp_out/"
        notify("wgs_snp 全部完成 ✅",
               f"台灣吳郭魚耐鹽 WGS/SNP 分析完成！\n\n{stats}",
               tags="tada")
        print("[monitor] DONE", flush=True)
        break

    # Error detection (throttle to once per 5 min)
    if ERROR_RE.search(line) and (time.time() - last_error) > 300:
        last_error = time.time()
        notify("wgs_snp 發生錯誤 ❌",
               f"偵測到錯誤：\n{line}\n\nLog: {LOG_FILE}",
               tags="warning")
        print(f"[notify] ERROR: {line}", flush=True)
