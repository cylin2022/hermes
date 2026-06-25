# Hermes 流程開發前置確認清單

每次開發新 Snakemake workflow 或修改現有 rule 時，依序完成以下各項確認再進行完整執行。

---

## 核心原則：「能跑完」≠「輸出正確」

> **Pilot test 驗收的是執行路徑，不是科學品質。**

這個清單通過了，只代表程式結構沒有崩潰。實際上有三類問題是測試無法自動抓到的：

| 問題類型 | 為何測試沒抓到 | 解法 |
|---------|-------------|------|
| **中間產物錯誤**（header 數量、BAM 品質、重複率偏差） | Pilot test 只確認「有輸出」，沒確認「輸出值合理」 | 第五關加**品質數字表**，必須人眼確認 |
| **Shell/Python 互動 bug**（`\1` → SOH、sed 靜默錯誤） | 程式不報錯，錯誤藏在輸出內容 | 第二關加**shell 命令預先獨立測試** |
| **操作程序問題**（中途改 Snakefile） | 測試環境沒有「production 期間有人在操作」這個維度 | 技術防護：啟動時凍結 Snakefile 副本（見第五關） |

**實際案例（2026-06）**：
- genome_annotation：`\1` 在 Snakemake shell block 中被 Python 解釋為 SOH（`\x01`），sed 靜默產出 1 個 header，39 個變成 1 個。Pilot test 沒有 `grep -c "^>"` 驗收。
- wgs_snp：`fixmate` 排序 bug 在 2 樣本 pilot test 沒有查 flagstat，流入 60 樣本後才被發現，觸發 ~11 天重跑。

**原則**：每個 workflow 的 pilot test 結束後，必須輸出一張**品質數字表**，讓人眼確認數字是否在合理範圍。通過數字驗收才算真正 pass。

---

## 第零關：新機器安裝時一次性設定

每台新主機安裝完 Hermes 後執行一次，生成針對該機器優化的 config：

```bash
# 偵測硬體，輸出建議參數
python3 hermes_configure.py

# 直接寫出設定檔
python3 hermes_configure.py --workflow wgs_snp --output my_config.yaml
```

檢查項目：

```
[ ] CPU cores / RAM / GPU 偵測結果合理
[ ] Docker 已安裝，且目前使用者在 docker group
[ ] NVIDIA Container Toolkit 可正常傳遞 GPU 進容器
[ ] 生成的 threads / dv_shards / mem_gb 已寫入實際使用的 config.yaml
```

**GPU job 排程原則**：
- CPU 工具（BWA-MEM2, samtools）用 `threads: N` 控制
- GPU 工具（DeepVariant）用 `threads: 1` + `resources: gpu=1`
- Snakemake 啟動時加 `--resources gpu=N`（N = GPU 數量）
- GPU job 的並行度完全由 GPU 資源決定，不受 CPU slots 影響

---

## 第一點五關：RNA 證據確認（genome_annotation 啟動前必查）

**BRAKER3 在沒有 RNA 證據的情況下，基因模型品質顯著下降**：外顯子邊界依賴遠緣物種蛋白質推斷，偽陽性基因數上升。

啟動 genome_annotation workflow 前，先確認：

```
[ ] 是否有 IsoSeq HiFi 資料？
      → 有：config 設定 isoseq_dir 指向 *.hifi_reads.fastq.gz 所在目錄
      → 無：是否有短讀長 RNA-seq 已映射的 BAM？→ 設定 rnaseq_bam
      → 都沒有：在 config 備注中記錄「蛋白質僅供第一版草稿」

[ ] isoseq_threads × 樣本數 ≤ 可用 CPU cores
      (8 樣本 × 16 threads = 128 cores；本機 160 cores OK)

[ ] isoseq_dir 內的 fastq.gz 命名格式為 {sample_name}.hifi_reads.fastq.gz

[ ] protein_fasta 已設定（Metazoa.fa 或物種相近的蛋白質集合）
```

**教訓（TWN_Hamaguri 2026-06-17）**：IsoSeq 資料已存在但尚未映射，為趕進度以蛋白質僅版啟動 BRAKER3。事後需重跑 BRAKER3 + DIAMOND + InterProScan，浪費 ~4 天運算。等映射完成（~3 h）遠比重跑划算。

---

## 第一關：Docker 容器規則（每個 docker run 都要確認）

```
[ ] --user $(id -u):$(id -g)  已加入每個 docker run
[ ] 輸出目錄的 -v bind-mount 路徑已用 realpath 解析（避免符號連結問題）
[ ] 確認 output 目錄在 container 啟動前已存在（mkdir -p）
```

**背景：** Docker 預設以 root 執行。若輸出檔案為 root 所有，Snakemake 的 `check_and_touch_output()` 在呼叫 `os.utime()` 時會 PermissionError，導致 job 標記失敗甚至 crash，且輸出檔案無法被 tabix 等後續工具覆寫。

---

## 第二關：工具行為確認（整合入 pipeline 之前）

針對每個新工具，手動執行一次單步驟測試：

```bash
# 取最小染色體做測試（快速重現）
CHROM="NC_XXXXXX.1"   # 換成你的最小 scaffold/chr

# 範例：DeepVariant 單獨測試
samtools view -b full.bam $CHROM > test.bam
samtools index test.bam
docker run --rm --gpus all --user $(id -u):$(id -g) \
    -v $(dirname $(realpath test.bam)):$(dirname $(realpath test.bam)) \
    -v $(dirname $(realpath genome.fa)):$(dirname $(realpath genome.fa)):ro \
    -v $(realpath outdir):$(realpath outdir) \
    google/deepvariant:1.10.0-gpu \
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=genome.fa \
        --reads=test.bam \
        --output_gvcf=outdir/test.g.vcf.gz \
        --regions $CHROM \
        --num_shards=4
```

確認項目：

```
[ ] exit code 是否如預期（工具說明書中的正常退出碼）
[ ] 輸出檔案的 owner 是 cylin（不是 root）
[ ] 輸出格式符合下游 rule 的預期（壓縮格式、index 有無）
[ ] 工具是否有隱性的輸入排序需求（fixmate 需 name-sort；MarkDuplicates 需 coord-sort）
[ ] 工具是否需要輔助文件（SnpEff: cds.fa/protein.fa；STAR: genome index）
[ ] 用 find 確認輸出檔案的「實際路徑」再寫入 Snakemake output:（很多工具會在子目錄輸出）
      find {outdir} -name "*.TElib.fa" -o -name "*.gtf" | head -5
```

**shell 命令預先獨立測試（加入 Snakemake 前必做）：**

Snakemake shell block 是 Python 字串，有非直覺的轉義規則。任何含有 `\`、`$`、`{}` 的 sed/awk 命令，**先在 bash 獨立測試一行**，確認輸出，再放進 Snakemake：

```bash
# 範例：先在 bash 確認 sed back-reference 效果
echo ">chr1 scaffold length=100" | sed 's/^\(>[^ ]*\) .*/\1/'
# 確認輸出是 ">chr1" 之後，Snakemake shell block 中改為 \\1

# 範例：確認 awk 變數傳遞
echo "test" | awk -v OFS='\t' '{print $1, NR}'
```

核心規則：
- `\1` → Snakemake shell block 中必須寫 `\\1`（Python 解釋 `\\1` 為字面 `\1`，sed 再解釋為 back-reference）
- `\\t` → shell block 中寫 `\\\\t`（兩層轉義）
- 靜默錯誤的症狀：輸出比預期少、行數減少、不可見字元開頭

**samtools pipeline 排序需求速查（違反者導致無聲錯誤）：**

| 工具 | 輸入需求 | 驗證指令 |
|------|---------|---------|
| `samtools fixmate -m` | **name-sorted** | `samtools view -H bam \| grep SO:` 應為 `queryname` |
| `samtools markdup` | **name-sorted**（fixmate 後接） | 同上 |
| `samtools sort`（最終） | 無限制（輸出 coord-sorted） | — |
| `samtools index` | **coord-sorted** | `samtools quickcheck bam` |
| DeepVariant | **coord-sorted + indexed** | `samtools quickcheck bam` |
| GATK MarkDuplicates | **coord-sorted** | `samtools view -H bam \| grep SO:` 應為 `coordinate` |

**教訓（wgs_snp 2026-06 事故）**：BWA-MEM2 多執行緒輸出未必 name-sorted，直接接 `fixmate` 在大樣本下會有亂序 read pair，導致 mate score 不準確、markdup 重複識別偏差。正確順序：`bwa-mem2 | samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup`。此 bug 在小測試資料（少 reads、低 thread）時不會觸發，難以察覺。

---

## 第三關：Snakemake rule 設計規範

```
[ ] shell block 開頭加 set -euo pipefail
[ ] 每個 tool 步驟後加明確輸出驗證：
        test -s {output.xxx} || { echo "ERROR: empty output" >&2; exit 1; }
[ ] || true 只用在「已確認工具正常完成但 exit 非零」的情況，並加上 comment 說明原因
[ ] 宣告 output: 的所有檔案名稱與 shell 實際產生的一致（包含副檔名大小寫）
[ ] 工具如有 temp/intermediate 目錄，在 shell 末尾 rm -rf 清除
```

---

## 第四關：設定檔與模板同步

每次修改工具版本或新增參數後，同步更新以下三個地方：

| 更動內容 | 需同步更新 |
|---------|-----------|
| Docker image 版本號 | `Snakefile`、`config_template.yaml`、`INSTALL.md` |
| 新工具/新步驟 | `config_template.yaml`（新增參數）、`README.md`（更新工具表） |
| 輸出副檔名改變 | 所有引用該輸出的下游 rule input |

確認指令：

```bash
# 確認 config_template.yaml 中的版本號與 Snakefile 一致
grep "dv_docker_image" workflows/wgs_snp/config_template.yaml workflows/wgs_snp/Snakefile
```

---

## 第五關：小型前導測試（前兩個樣本 × 最小染色體）

在正式跑全部樣本前，先執行小型測試：

### 5.1 建立測試用 samplesheet

```bash
# 取前 2 個樣本
head -3 samplesheet.csv > test_samplesheet.csv   # header + 2 samples
```

### 5.2 建立測試用 genome（只保留最小的 3 條 scaffold）

```bash
# 找最小的 3 條 scaffold
samtools faidx genome.fa
sort -k2,2n genome.fa.fai | head -3 | awk '{print $1}' > test_chroms.txt
samtools faidx genome.fa $(cat test_chroms.txt | tr '\n' ' ') > test_genome.fa
samtools faidx test_genome.fa
```

### 5.3 修改 config 指向測試資料

```yaml
samplesheet: "test_samplesheet.csv"
genome_fasta: "test_genome.fa"
outdir: "test_run_output"
```

### 5.4 執行前導測試

```bash
# 先 dry-run 確認 DAG
snakemake -s Snakefile --cores 8 --use-conda \
    --configfile test_config.yaml --dry-run

# 確認無誤後執行（不超過 8 cores，避免佔用生產資源）
snakemake -s Snakefile --cores 8 --use-conda \
    --configfile test_config.yaml
```

### 5.5 前導測試驗收

```
[ ] 所有 rule 完成，無 error（包含 SnpEff annotation rule）
[ ] 輸出檔案 owner 均為 cylin
[ ] 輸出 VCF/BAM 可被 tabix/samtools 正常讀取
[ ] Log 無 PermissionError / OSError / utime 相關訊息
[ ] Snakemake 末尾顯示 "N of N steps (100%) done"
```

**WGS 流程額外驗收（BAM 品質指標，不可省略）：**

```bash
# 1. 確認 BAM 排序正確
samtools view -H test_out/alignment/*.markdup.bam | grep "^@HD"
# 應顯示 SO:coordinate（最終 markdup BAM）

# 2. 確認 flagstat 重複率合理
samtools flagstat test_out/alignment/*.markdup.bam
# WGS 正常範圍：重複率 5–25%
# 若 >40% → markdup 輸入排序錯誤（name-sort 問題）
# 若 <1%  → markdup 可能未正確執行

# 3. 確認 fixmate 有效（mate score 欄位存在）
samtools view test_out/alignment/*.namesorted.bam | head -5 | cut -f12-
# 應有 ms:i: 欄位（mate score）；若無 → fixmate 輸入未 name-sorted
```

**教訓（wgs_snp 2026-06 事故）**：未在 pilot test 驗收 flagstat 重複率，導致 `fixmate` 排序 bug 流入全量 60 樣本執行。Bug 修正後全部 60 個 BAM 被重建（5 天 CPU），49 個樣本觸發 DeepVariant 重跑（6 天 GPU）。若 pilot test 就確認重複率，代價僅為 2 個樣本的重跑。

### 5.6 中間產物品質數字確認（各 workflow 必填表）

Pilot test 完成後，填寫以下表格並**人眼確認數字合理**。不合理的數字必須追查原因，不能直接進入全量執行。

**wgs_snp：**

```bash
SAMPLE="test_sample"   # 替換成實際測試樣本名稱

echo "=== WGS SNP Pilot 品質確認 ==="
echo "BAM 排序："
samtools view -H test_run_output/alignment/${SAMPLE}.markdup.bam | grep "^@HD"
# 預期：SO:coordinate

echo "重複率："
samtools flagstat test_run_output/alignment/${SAMPLE}.markdup.bam | grep -E "duplicate|total"
# 預期：重複率 5–25%；遠低於 5% = markdup 未執行；遠高於 25% = 排序錯誤

echo "mate score 欄位："
samtools view test_run_output/alignment/${SAMPLE}.markdup.bam | head -5 | cut -f12- | grep -o "ms:i:[0-9]*" | head -3
# 預期：有 ms:i: 欄位；若無 = fixmate 輸入未 name-sorted

echo "gVCF 記錄數："
bcftools stats test_run_output/gvcf/${SAMPLE}.g.vcf.gz | grep "^SN" | head -5
# 預期：有合理數量的 records（小染色體 ~數千至數萬）
```

**genome_annotation：**

```bash
echo "=== Genome Annotation Pilot 品質確認 ==="
echo "FASTA header 數："
grep -c "^>" test_run_output/genome_shortid.fa
# 預期：與原始 genome.fa header 數一致

echo "Repeat library 大小："
grep -c "^>" test_run_output/repeat_library.fa
# 預期：>100 個 repeat family（小基因組可接受 50+）

echo "BRAKER3 gene 數量（若已完成）："
grep -c "transcript_id" test_run_output/annotation/braker.gtf
# 預期：蛋白質編碼基因 >5000（脊椎動物），軟體動物 >10000
```

**rnaseq / scrnaseq：**

```bash
echo "=== RNA-seq Pilot 品質確認 ==="
echo "STAR mapping rate："
grep "Uniquely mapped reads %" test_run_output/star_log/${SAMPLE}_Log.final.out
# 預期：>60%（高品質）；<40% 需查基因組版本或 reads 品質

echo "featureCounts 分配率："
grep "Successfully assigned" test_run_output/counts/${SAMPLE}.summary
# 預期：>50%；若低 = strand 設定錯誤或 annotation 不匹配
```

| 驗收項目 | 預期範圍 | 實際值 | Pass? |
|---------|---------|--------|-------|
| BAM SO: tag | coordinate | | |
| 重複率 | 5–25% | | |
| ms:i: 欄位存在 | 是 | | |
| FASTA header 數 | = 原始數量 | | |
| mapping rate | >60% | | |

**所有欄位 Pass 才算 pilot test 通過。**

---

## 第六關：監控設定確認

```
[ ] snakemake.pid 在啟動後存在且 PID 確實是活著的程序
    ps -p $(cat runs/<run_id>/snakemake.pid)
[ ] notify.py / hermes_cmd.py 的 LOG_FILE 路徑指向正確的 run
[ ] 通知測試：確認 ntfy / email 能收到
[ ] mail_watcher 已在背景執行（pgrep -f mail_watcher 或 monitor_wgs.py）
```

---

## 常見工具已知問題速查

| 工具 | 問題 | 解決方式 |
|-----|-----|---------|
| DeepVariant ≤ 1.6.1 | exit 1 even on success | 升級至 1.10.0-gpu 或更新版 |
| DeepVariant 1.10.0-gpu | 內部 tabix 失敗 exit 1 | `\|\| true`；外部重建 tabix index |
| SnpEff + NCBI GTF | 缺 cds.fa/protein.fa → build exit 255 | 加 `-noCheckCds -noCheckProtein` |
| SnpEff | 只認 `genes.gff`，不認 `genes.gff3` | 副檔名用 `gff`（不是 `gff3`） |
| **samtools fixmate -m（重要）** | 需要 name-sorted input；BWA-MEM2 多執行緒輸出不保證 name-sorted → fixmate mate score 錯誤 → markdup 重複識別偏差。小測試不會觸發（reads 少時 BWA 輸出幾乎有序），大樣本才暴露。**驗證**：`samtools view -H bam \| grep "SO:"` 應為 `queryname`；`samtools flagstat` 重複率應 5–25% | 正確順序：`bwa-mem2 \| samtools sort -n \| samtools fixmate -m \| samtools sort \| samtools markdup`；pilot test 必查 flagstat 重複率 |
| **Snakemake 中途 code change 連鎖重跑** | pipeline 執行中修改 Snakefile rule → Snakemake 預設 code-change 觸發 → 受影響 rule 全樣本重跑 → 下游時間戳更新 → 進一步觸發後續 rule 重跑（wgs_snp 2026-06：60 樣本 BAM 全重建，49 樣本 DeepVariant 重跑，損失 ~11 天 GPU+CPU） | 修改前先 dry-run 評估影響範圍；style 改動用 `--rerun-triggers mtime` 重啟；bug fix 改動評估代價後再決定；生產期儘量不改 Snakefile |
| Docker 任意工具 | 輸出為 root-owned | 加 `--user $(id -u):$(id -g)` |
| DeepVariant _tmp 殘留 | 前次沒有 `--user` 的 tmp 中有 root-owned tfrecord；新 DV 以 cylin 跑時無法覆寫 → g.vcf.gz 不產出 | shell block 開頭加 `rm -rf {params.tmp_dir}` 再 `mkdir -p` |
| Snakemake crash 後重啟 | 目錄鎖定 | `snakemake --unlock` 再重啟 |
| Snakemake 在 conda 外執行 | command not found | `source ~/miniforge3/etc/profile.d/conda.sh && conda activate base` |
| RepeatModeler 中斷 | 重跑從頭 | `-recoverDir RM_*/` 接續 |
| **RepeatModeler** | **基因組 >1.5 GB 時呈指數成長**：2 GB 基因組 Round 5 = 22.8M pairwise BLAST，Round 6 ≈ 200M，可能跑數週 | 改用 **EDTA**（結構偵測，無 all-vs-all BLAST，~24-48 h）；設 `repeat_tool: "edta"` |
| **BRAKER3 conda** | `braker3=3.0.8` conda env 缺 `Scalar::Util::Numeric` Perl 模組 → `Can't locate Scalar/Util/Numeric.pm` exit 255 | 改用 Docker：`teambraker/braker3:latest`（已確認有該模組） |
| **braker4-allinone:latest** | 即使以 root 執行也缺 `Scalar::Util::Numeric`（conda Perl 5.32 內建版本不含此模組） | **不要用此 image**；改用 `teambraker/braker3:latest` |
| **BRAKER3 Docker + `--user`** | `--user $(id -u):$(id -g)` 導致容器內 Perl 模組路徑 permission denied | BRAKER3 必須以 **root** 在容器內執行，**不加 `--user`** |
| **BRAKER3 workingdir 權限** | Container 以 root 寫入，但 workingdir 預設 755，BRAKER3 回報 `Do not have write permission for /workingdir` | docker run 前先 `chmod 777 {params.outdir}` |
| **BRAKER3 輸出檔名** | 有時輸出 `augustus.hints.aa` 而非 `braker.aa`，導致 Snakemake missing output error | shell block 末尾加 fallback：`if [ ! -f braker.aa ] && [ -f augustus.hints.aa ]; then cp augustus.hints.aa braker.aa; fi` |
| **BRAKER3 root-owned 輸出** | Docker 以 root 執行，輸出目錄所有檔案為 root-owned；`sudo chown` 需要密碼靜默失敗 | 用 busybox 容器修正：`docker run --rm -v <outdir>:/workingdir busybox sh -c "chown -R uid:gid /workingdir"`；**必須在 Snakemake shell block 內執行**（不是事後手動），否則 Snakemake postprocess `os.utime()` 會拋 `PermissionError: [Errno 13]` 導致 job 標記失敗 |
| **BRAKER3 root-owned → Snakemake crash 後復原** | braker.gtf/braker.aa 已產出但仍為 root-owned，Snakemake 標記 job 失敗無法繼續 | 手動執行 busybox chown 修正擁有者 → `--rerun-incomplete` 重啟；無需重跑 BRAKER3 |
| **EDTA output 路徑錯誤** | EDTA 實際輸出在 `genome.mod.EDTA.final/genome.mod.EDTA.TElib.fa`（有一層子目錄），Snakefile 若只寫 `genome.mod.EDTA.TElib.fa` → 永遠找不到 output，不斷重跑 | Snakemake `output:` 必須包含 `.mod.EDTA.final/` 子目錄層；修正 Snakefile 後若輸出已存在，用 `--rerun-triggers mtime` 重啟避免因 code-change 觸發重跑 |
| **BUSCO v5 vs v6 / odb10 vs odb12** | `busco=5.7.1` 不支援 odb12 資料集 → `ERROR: BUSCO v5 only works with datasets from OrthoDB v10` | `envs/busco.yaml` 改為 `busco=6.0.0`；odb12 lineage 必須搭配 BUSCO v6 |
| **Snakemake `if X:` 空 body** | `if condition:` 後只有 comment 沒有程式碼 → `IndentationError: expected an indented block` | 合併條件：`if A and B:` / `if A and not B:`，避免空 if body |
| **nohup 腳本 conda not found** | nohup 啟動時 PATH 不含 miniforge3 → `conda: command not found` | nohup 腳本開頭加 `export PATH="/home/cylin/miniforge3/bin:$PATH"` |
| **Snakemake shell 反斜線轉義（`\1` → `\x01`）** | Snakemake shell block 是 Python 字符串：`\1` 被 Python 解釋為 SOH 控制字元（`\x01`），實際傳給 shell 的是 SOH 而非 sed 的 back-reference `\1`。症狀：log 顯示替換結果為空，或輸出檔案以不可見控制字元開頭，`grep "^>"` 結果大幅減少（TWN_Hamaguri 案例：39 個 FASTA header 變成 1 個）。 | 在 Snakemake shell block 中永遠用 `\\1`（Python 解釋為字面 `\1`，再傳給 sed 作 back-reference）。驗證方式：啟動後立即 `grep -c "^>" output.fa` 確認 header 數量與 input 一致。 |

---

## 中途修改 Snakefile 的風險與處置

**核心問題**：Snakemake 預設使用 `code` 觸發器——只要 rule 的 shell block 或參數有任何改動，所有依賴該 rule 的下游樣本都會被標記為需重跑。若在全量執行中途修改 Snakefile，可能觸發數十個樣本的連鎖重跑，浪費數天運算。

**教訓（wgs_snp 2026-06 事故）**：

```
bwa_mem2 rule 修正（加 samtools sort -n）
  → Snakemake code-change 觸發全部 60 樣本重跑 bwa_mem2
    → 新 BAM 時間戳比舊 gVCF 新
      → 49 個樣本排隊重跑 DeepVariant（1 GPU，耗時 ~6 天）
```

### 修改前必問三個問題

```
1. 這個修改是 bug fix（影響結果）還是 code style（不影響結果）？
   → bug fix：接受重跑，但評估影響範圍再決定時機
   → code style：用 --rerun-triggers mtime 重啟，跳過 code-change 觸發

2. 目前執行到哪一步？已完成的步驟有多少？
   → 若大部分樣本都已過了受影響的 rule → 等全數完成後再修改
   → 若只有少數完成 → 修改後重跑代價較小，可接受

3. 修改後重跑需要多少時間？是否有資源支撐？
   → 估算：受影響樣本數 × 單樣本耗時 × 資源瓶頸（GPU/CPU）
   → 若超過 24 小時 → 一定要先評估，再決定是否修改
```

### 依情境選擇重啟方式

| 情境 | 重啟指令 | 說明 |
|------|---------|------|
| bug fix，需重跑 | `--rerun-incomplete` | 重跑有 code change 或未完成的 job |
| code style，不需重跑 | `--rerun-triggers mtime` | 只有 input 比 output 新才重跑，忽略 code change |
| 輸出已正確但時間戳過舊 | `touch` 輸出檔案後 `--rerun-triggers mtime` | 告訴 Snakemake 輸出是最新的 |

```bash
# code style 修改後安全重啟（不觸發 code-change 重跑）
nohup snakemake -s Snakefile \
    --configfile config.yaml \
    --cores 128 --use-conda \
    --rerun-incomplete \
    --rerun-triggers mtime \   # ← 關鍵：只看時間戳，忽略 code change
    > snakemake.log 2>&1 &

# 確認哪些 job 會被觸發（先 dry-run）
snakemake --dry-run --rerun-triggers mtime ... 2>&1 | grep "^Job\|^rule" | head -30
```

### 預防：開發與生產分離

- **開發期**（pipeline 設計階段）：修改 Snakefile → pilot test → 確認無誤
- **生產期**（全量執行中）：除非有嚴重 bug，否則**不修改正在執行的 Snakefile**
- 若確實需要修改：先等當前批次跑完，再修改、再以正確 flag 重啟

---

## 崩潰後標準恢復程序

```bash
# 1. 確認 Snakemake 確實已停止
ps aux | grep snakemake

# 2. 解除鎖定
snakemake -s workflows/<wf>/Snakefile --unlock

# 3. 確認問題已修復後，以 --rerun-incomplete 重啟
source ~/miniforge3/etc/profile.d/conda.sh && conda activate base
nohup snakemake -s workflows/<wf>/Snakefile \
    --cores 128 --use-conda --rerun-incomplete \
    --configfile <config.yaml> \
    > runs/<run_id>/snakemake.log 2>&1 &

# 4. 記錄新 PID
echo $! > runs/<run_id>/snakemake.pid
echo $! > runs/<run_id>/pid
```
