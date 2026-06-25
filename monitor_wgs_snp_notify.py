#!/usr/bin/env python3
"""
監控 WGS SNP pipeline 完成，整理結果並 email 通知。
等待策略：Snakemake PID 結束（確保 T-13-1 重跑後的最終結果）。
"""
import subprocess
import sys
import time
import os
from pathlib import Path

SNAKEMAKE_PID = 1152019
OUTDIR    = Path("/home/cylin/WGS_AntiSalt_tilapia/wgs_snp_out")
ANNOTATED = OUTDIR / "annotation" / "snps.annotated.vcf.gz"
SUMMARY   = OUTDIR / "stats" / "variant_summary.txt"
SNP_VCF   = OUTDIR / "vcf" / "filtered" / "snps.PASS.vcf.gz"
INDEL_VCF = OUTDIR / "vcf" / "filtered" / "indels.PASS.vcf.gz"
HERMES    = Path("/home/cylin/hermes")

sys.path.insert(0, str(HERMES))


def sh(cmd: str) -> str:
    try:
        return subprocess.check_output(cmd, shell=True, text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "?"


def fmt(n: str) -> str:
    try:
        return f"{int(n):,}"
    except Exception:
        return n


def mb(p: Path) -> str:
    return f"{p.stat().st_size / 1e6:.1f} MB" if p.exists() else "N/A"


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def wait_for_snakemake(pid: int, interval: int = 120):
    print(f"[wgs_notify] 等待 Snakemake PID {pid} 結束...", flush=True)
    while pid_alive(pid):
        time.sleep(interval)
    print(f"[wgs_notify] Snakemake PID {pid} 已結束", flush=True)


def parse_summary() -> dict:
    if not SUMMARY.exists():
        return {}
    return {
        "n_samples":   sh(f"grep 'number of samples' {SUMMARY} | head -1 | awk '{{print $NF}}'"),
        "snp_count":   sh(f"awk '/PASS SNPs/{{f=1}} f && /number of records/{{print $NF; exit}}' {SUMMARY}"),
        "indel_count": sh(f"awk '/PASS INDELs/{{f=1}} f && /number of records/{{print $NF; exit}}' {SUMMARY}"),
        "tstv":        sh(f"grep '^TSTV' {SUMMARY} | head -1 | awk '{{print $5}}'"),
        "multi":       sh(f"awk '/PASS SNPs/{{f=1}} f && /multiallelic sites/{{print $NF; exit}}' {SUMMARY}"),
    }


def snpeff_impacts() -> dict:
    if not ANNOTATED.exists():
        return {}
    print("[wgs_notify] 統計 SNPeff impact（1–2 分鐘）...", flush=True)
    high     = sh(f"bcftools view -H {ANNOTATED} | cut -f8 | grep -c '|HIGH|' || echo 0")
    moderate = sh(f"bcftools view -H {ANNOTATED} | cut -f8 | grep -c '|MODERATE|' || echo 0")
    return {"high": high, "moderate": moderate}


def main():
    # 等 Snakemake 結束（T-13-1 → GLnexus → filter → SNPeff 全部跑完）
    wait_for_snakemake(SNAKEMAKE_PID)

    print("[wgs_notify] 收集統計...", flush=True)
    s   = parse_summary()
    imp = snpeff_impacts()

    lines = [
        "【鹽度 WGS SNP pipeline 完成】",
        "",
        "▌ 基本統計",
        f"  樣本數       : {s.get('n_samples', '?')}",
        f"  PASS SNPs    : {fmt(s.get('snp_count', '?'))}",
        f"  PASS INDELs  : {fmt(s.get('indel_count', '?'))}",
        f"  Ts/Tv ratio  : {s.get('tstv', '?')}  （WGS 品質基準 >1.8）",
        f"  多等位位點   : {fmt(s.get('multi', '?'))}",
        "",
    ]

    if imp:
        lines += [
            "▌ SNPeff 功能影響分類",
            f"  HIGH impact     : {fmt(imp.get('high', '?'))}  （停止密碼子、移碼突變等）",
            f"  MODERATE impact : {fmt(imp.get('moderate', '?'))}  （胺基酸替換等）",
            "",
        ]

    lines += [
        "▌ 輸出檔案",
        f"  SNPs PASS    : {SNP_VCF}",
        f"               ({mb(SNP_VCF)})",
        f"  INDELs PASS  : {INDEL_VCF}",
        f"               ({mb(INDEL_VCF)})",
    ]
    if ANNOTATED.exists():
        lines += [
            f"  SNPeff 注釋  : {ANNOTATED}",
            f"               ({mb(ANNOTATED)})",
        ]
    lines += [
        f"  統計摘要     : {SUMMARY}",
        "",
        "▌ 建議後續分析",
        "  1. snp_association workflow（GWAS / GEMMA LMM / Fst 掃描）",
        "  2. PLINK2 PCA — 確認耐鹽 vs 非耐鹽群體結構",
        "  3. HIGH impact SNPs 候選基因清單",
        "",
        "— 小賀 自動通知",
    ]

    body    = "\n".join(lines)
    subject = f"WGS SNP 完成：{fmt(s.get('snp_count','?'))} SNPs，{s.get('n_samples','?')} 樣本"

    import notify
    notify.send(subject, body)
    print(f"[wgs_notify] Email 已寄出：{subject}", flush=True)
    print(body, flush=True)


if __name__ == "__main__":
    main()
