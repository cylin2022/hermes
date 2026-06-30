# diagnose-snakemake-failure

> **Skill**: 系統化診斷 Hermes Snakemake workflow 失敗
> **觸發時機**: 任何 workflow 失敗、rule error、crash 後恢復

---

## 第一步：取得錯誤全貌

```
get_log(run_id, lines=100)
```

掃描關鍵訊號（按優先順序）：

| 關鍵字 | 錯誤類型 | 跳至 |
|--------|---------|------|
| `Error: directory is locked` | Snakemake 未正常結束 | § Snakemake 鎖定 |
| `PermissionError: [Errno 13]` / `os.utime` | Docker root-owned 輸出 | § Docker 權限錯誤 |
| `MemoryError` / `Killed` / `OOM` | 記憶體不足 | § OOM |
| `MissingOutputException` | 輸出檔案未產出 | § 缺少輸出 |
| `conda: command not found` | conda 環境未載入 | § conda 環境問題 |
| `IndentationError` / `SyntaxError` | Snakefile 語法錯誤 | § 語法錯誤 |
| `code-change` / rule 被重新觸發 | 中途改 Snakefile | § Code-change 連鎖 |
| `Can't locate Scalar/Util/Numeric.pm` | BRAKER3 Perl 模組缺失 | § BRAKER3 特定 |
| `exit 1` (DeepVariant) | DV 已知問題 | § DeepVariant 特定 |
| header 數量異常減少 | `\1` 轉義 bug | § Shell 轉義 |

---

## § Snakemake 鎖定

**症狀**：`Error: directory is locked`

```bash
# 確認 Snakemake 確實已停止
ps aux | grep snakemake

# 解鎖（只在確認 Snakemake 不在執行時才操作）
snakemake -s workflows/<wf>/Snakefile --unlock

# 重啟
source ~/miniforge3/etc/profile.d/conda.sh && conda activate base
nohup snakemake -s workflows/<wf>/Snakefile \
    --cores 128 --use-conda --rerun-incomplete \
    --configfile <config.yaml> \
    > runs/<run_id>/snakemake.log 2>&1 &
echo $! > runs/<run_id>/snakemake.pid
```

---

## § Docker 權限錯誤

**症狀**：`PermissionError: [Errno 13]` 或輸出檔案為 root 所有

**根本原因**：Docker 未加 `--user $(id -u):$(id -g)`，輸出為 root-owned，
Snakemake 呼叫 `os.utime()` 時 PermissionError。

**例外**：BRAKER3 必須以 root 在容器內執行（不加 `--user`）。
BRAKER3 root-owned 輸出需用 busybox 容器修正：

```bash
# BRAKER3 輸出修正（必須在 Snakemake shell block 內，不能事後手動）
docker run --rm \
  -v <outdir>:/workingdir \
  busybox sh -c "chown -R $(id -u):$(id -g) /workingdir"
```

**其他工具**：在 `docker run` 加入 `--user $(id -u):$(id -g)`

**已損壞狀態恢復**：

```bash
# 手動修正已存在的 root-owned 輸出
docker run --rm -v <outdir>:/workingdir busybox \
  sh -c "chown -R <uid>:<gid> /workingdir"
# 然後 --rerun-incomplete（無需重跑 BRAKER3）
```

---

## § OOM（記憶體不足）

**症狀**：`MemoryError` / `Killed` / dmesg 中 OOM killer

**策略**（依序嘗試）：

1. 確認當前記憶體使用：
   ```bash
   free -h
   ```

2. 減少 config 中的 `mem_gb`（初始建議減半）

3. 減少 `threads`（少核心 = 少並行記憶體占用）

4. DeepVariant 特別做法：減少 `dv_shards`
   ```yaml
   dv_shards: 16   # 從 64 降到 16
   ```

5. GLnexus：增加 `--mem-gbytes` 限制或分批執行

> 注意：本伺服器有 2.2 TiB RAM，真正的 OOM 通常是
> 單 rule 的 memlimit 設定過低，而非全系統記憶體不足。
> 先查 rule 的 `resources: mem_mb` 設定。

---

## § 缺少輸出（MissingOutputException）

**診斷步驟**：

```bash
# 1. 確認檔案是否真的不存在
ls -la <expected_output_path>

# 2. 確認 output 目錄是否存在（很多工具需要預先 mkdir）
ls <output_dir>

# 3. 確認 Docker 容器輸出路徑（bind-mount 問題）
# 工具輸出路徑可能在子目錄，例如：
# EDTA 輸出在 genome.mod.EDTA.final/genome.mod.EDTA.TElib.fa
# 而非根目錄的 genome.mod.EDTA.TElib.fa

# 4. 確認工具 exit code（log 末尾）
tail -50 runs/<run_id>/snakemake.log | grep -E "exit|Error|error"
```

**常見原因**：
- EDTA：輸出在子目錄 `.mod.EDTA.final/`，Snakefile output: 路徑漏掉這層
- BRAKER3：有時輸出 `augustus.hints.aa` 而非 `braker.aa`
  → shell block 末尾需有 fallback：
  ```bash
  if [ ! -f braker.aa ] && [ -f augustus.hints.aa ]; then
      cp augustus.hints.aa braker.aa
  fi
  ```

---

## § conda 環境問題

**症狀**：`conda: command not found` 或 `command not found`

```bash
# nohup 腳本或新 shell 中需要先載入 conda
export PATH="/home/cylin/miniforge3/bin:$PATH"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate base
```

---

## § Code-change 連鎖重跑

**症狀**：修改 Snakefile 後，大量已完成的 job 被重新觸發

**評估影響**（先 dry-run）：
```bash
snakemake --dry-run --rerun-triggers mtime \
    -s workflows/<wf>/Snakefile \
    --configfile <config.yaml> 2>&1 | grep "^Job\|^rule" | head -30
```

**依情境選擇重啟方式**：

| 修改類型 | 重啟指令 |
|---------|---------|
| Bug fix（影響結果） | `--rerun-incomplete`（接受重跑） |
| Code style（不影響結果） | `--rerun-triggers mtime`（跳過 code-change） |
| 輸出正確但時間戳過舊 | `touch` 輸出檔 + `--rerun-triggers mtime` |

> [!WARNING]
> 教訓（wgs_snp 2026-06）：中途修改 `bwa_mem2` rule
> 觸發 60 樣本重跑 BAM、49 樣本重跑 DeepVariant，損失約 11 天。
> 生產執行中不要修改 Snakefile，除非有嚴重 bug。

---

## § BRAKER3 特定問題

**`Can't locate Scalar/Util/Numeric.pm`**

- 原因：`braker3` conda env 或 `braker4-allinone:latest` 缺少此 Perl 模組
- 解法：改用 `teambraker/braker3:latest` Docker image
- 不要用 `braker4-allinone:latest`（即使以 root 執行也缺此模組）

**BRAKER3 workingdir 權限**

```bash
# Docker run 前必須
chmod 777 <params.outdir>
```

**BRAKER3 root-owned 輸出 → Snakemake crash 後無法繼續**

```bash
docker run --rm -v <braker_outdir>:/workingdir busybox \
  sh -c "chown -R $(id -u):$(id -g) /workingdir"
# 修正後：--rerun-incomplete 重啟（不需重跑 BRAKER3）
```

---

## § DeepVariant 特定問題

**exit 1 即使成功**

- DV <= 1.6.1：exit 1 even on success → 升級至 1.10.0-gpu
- DV 1.10.0-gpu：內部 tabix 有時 exit 1
  → shell block 對 DV 步驟加 `|| true`，後面外部重建 tabix index：
  ```bash
  tabix -p vcf {output.gvcf}
  ```

**tmp 目錄殘留（root-owned）**

```bash
# shell block 開頭加
rm -rf {params.tmp_dir}
mkdir -p {params.tmp_dir}
```

---

## § Shell 轉義 Bug（`\1` → SOH）

**症狀**：FASTA header 數量大幅減少（例如 39 個變 1 個）

**原因**：Snakemake shell block 是 Python 字串，`\1` 被 Python 解釋為
SOH（`\x01`），而非 sed 的 back-reference `\1`

**修正**：在 Snakemake shell block 中永遠用 `\\1`

**驗證**：
```bash
grep -c "^>" output.fa
# 數量應與 input.fa 的 header 數相同
```

---

## § samtools 排序問題（WGS 流程）

**症狀**：markdup 重複率 > 40% 或 < 1%

**正確的 BAM 處理順序**：
```
bwa-mem2 | samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup
```

**驗證**：
```bash
# fixmate 前的 BAM 應為 name-sorted
samtools view -H *.namesorted.bam | grep "SO:"   # 應為 queryname

# 最終 markdup BAM 應為 coord-sorted
samtools view -H *.markdup.bam | grep "SO:"      # 應為 coordinate

# 重複率應在 5–25%
samtools flagstat *.markdup.bam | grep "duplicate"
```

---

## 標準恢復流程（任意 workflow）

```bash
# 1. 確認 Snakemake 已停止
ps aux | grep snakemake

# 2. 解鎖
snakemake -s workflows/<wf>/Snakefile --unlock

# 3. 重啟（問題已修復後）
source ~/miniforge3/etc/profile.d/conda.sh && conda activate base
nohup snakemake -s workflows/<wf>/Snakefile \
    --cores 128 --use-conda --rerun-incomplete \
    --configfile <config.yaml> \
    > runs/<run_id>/snakemake.log 2>&1 &
echo $! > runs/<run_id>/snakemake.pid
echo $! > runs/<run_id>/pid
```
