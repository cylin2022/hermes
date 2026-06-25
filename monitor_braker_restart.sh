#!/usr/bin/env bash
# 監控 BRAKER3 完成，自動重啟 Snakemake with --cores 160
# 目的：讓 BUSCO(40) + DIAMOND(80) + InterProScan(40) 同時跑，充分用滿 160 cores
set -euo pipefail

BRAKER_GTF="/home/cylin/TWN_Hamaguri/assembly_output/annotation/braker4/braker.gtf"
SNAKEMAKE_PID=1207445
CONFIGFILE="/home/cylin/hermes/runs/genome_annotation_hamaguri_edta_20260617/config.yaml"
SNAKEFILE="workflows/genome_annotation/Snakefile"
LOGDIR="/home/cylin/hermes/runs/genome_annotation_hamaguri_edta_20260617"
HERMES="/home/cylin/hermes"

export PATH="/home/cylin/miniforge3/bin:$PATH"

echo "[monitor_braker_restart] 啟動，等待 braker.gtf 出現..."
echo "[monitor_braker_restart] 監控: $BRAKER_GTF"

while true; do
    if [ -s "$BRAKER_GTF" ]; then
        echo "[monitor_braker_restart] braker.gtf 已就緒 ($(du -h $BRAKER_GTF | cut -f1))"

        # 等舊 Snakemake 自然退出（braker4 完成後它會繼續跑下游，但用的是舊設定）
        # 給它 60 秒處理 busybox chown 等收尾工作
        sleep 60

        # 停止舊 Snakemake（如果仍在執行）
        if kill -0 $SNAKEMAKE_PID 2>/dev/null; then
            echo "[monitor_braker_restart] 停止舊 Snakemake PID $SNAKEMAKE_PID..."
            kill $SNAKEMAKE_PID
            sleep 5
        fi

        # 解除 Snakemake 鎖定
        cd "$HERMES"
        python3 -m snakemake -s "$SNAKEFILE" \
            --configfile "$CONFIGFILE" \
            --unlock 2>/dev/null || true

        # 重啟：160 cores，讓 BUSCO+DIAMOND+InterProScan 並行
        RESTART_LOG="$LOGDIR/snakemake_postbraker_$(date +%Y%m%d_%H%M%S).log"
        nohup python3 -m snakemake \
            -s "$SNAKEFILE" \
            --configfile "$CONFIGFILE" \
            --cores 160 \
            --rerun-incomplete \
            --rerun-triggers mtime \
            --keep-going \
            -p \
            > "$RESTART_LOG" 2>&1 &
        NEW_PID=$!
        echo "[monitor_braker_restart] 新 Snakemake PID: $NEW_PID (--cores 160)"
        echo "[monitor_braker_restart] Log: $RESTART_LOG"
        echo "$NEW_PID" > "$LOGDIR/snakemake_postbraker.pid"

        # 通知
        python3 "$HERMES/notify.py" \
            "TWN_Hamaguri BRAKER3 完成，後處理啟動" \
            "braker.gtf 已就緒。重啟 Snakemake (PID $NEW_PID, --cores 160)。\nBUSCO(40) + DIAMOND(80) + InterProScan(40) 同時執行。" 2>/dev/null || true

        exit 0
    fi

    sleep 300  # 每 5 分鐘檢查一次
done
