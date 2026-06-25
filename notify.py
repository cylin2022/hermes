#!/usr/bin/env python3
"""Send email notification via Gmail SMTP. Called by monitor scripts."""
import json
import smtplib
import ssl
import sys
from email.mime.text import MIMEText
from pathlib import Path

CONFIG = Path.home() / ".hermes_notify.json"


def send(subject: str, body: str) -> None:
    cfg = json.loads(CONFIG.read_text())
    msg = MIMEText(body, "plain", "utf-8")
    msg["Subject"] = f"[小賀] {subject}"
    msg["From"]    = cfg["gmail_user"]
    msg["To"]      = cfg["notify_to"]

    ctx = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=ctx) as smtp:
        smtp.login(cfg["gmail_user"], cfg["gmail_app_password"])
        smtp.sendmail(cfg["gmail_user"], cfg["notify_to"], msg.as_string())


if __name__ == "__main__":
    subject = sys.argv[1] if len(sys.argv) > 1 else "通知"
    body    = sys.argv[2] if len(sys.argv) > 2 else ""
    send(subject, body)
    print(f"Sent: {subject}")
