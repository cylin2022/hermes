# intake-qc-checklist

> **Skill**: 新分析請求進來時，引導用戶確認資料品質，防止「跑完才發現樣本壞掉」
> **觸發時機**: 用戶描述新的分析需求、提供資料路徑、詢問「可以開始跑嗎？」

---

## 核心原則

> 「能跑完」≠「輸出正確」。Pilot test 通過只代表程式結構沒有崩潰。
> 在啟動任何 workflow 之前，先確認資料本身值得分析。

---

## 通用問題（所有 workflow 都問）

在規劃任何分析前，先確認：

```
□ 資料路徑確實存在？
  ls -la <r1_path> <r2_path>   # 確認檔案存在且非空

□ 檔案是否完整（非截斷）？
  gzip -t <file.fastq.gz>      # exit 0 = 完整；非零 = 損壞

□ 樣本名稱是否唯一、無空格、無特殊字元？
  # 空格和括號會讓 Snakemake wildcards 爆炸

□ 這批資料是否已跑過其他分析？結果在哪裡？
  # 避免重複計算，也確認資料「歷史」
```

---

## rnaseq（RNA-seq 差異表達）

### 資料進入前必問

```
1. 幾個樣本？幾個 condition？每個 condition 幾個 replicate？
   → 每組 < 3 replicate：DESeq2 統計力不足，需告知用戶結果信度低
   → 只有 1 replicate：無法做 DE，改跑 exploratory（PCA + heatmap）

2. 是否有成對設計（paired samples）？
   → 有：DESeq2 需加 design formula（~ batch + condition）

3. 資料來源（protocol）？
   → TruSeq Stranded（最常見）→ featureCounts -s 2
   → Unstranded → featureCounts -s 0
   → 不確定 → 啟動前先跑 infer_experiment.py（5 分鐘，值得）

4. 預估 library size 是否均一？
   → 請用戶提供 FastQC/multiqc 報告（若已有）
   → 若 library size 差異 > 10x → 建議先做 downsampling 再評估

5. 是否有 batch effect（不同批次、不同時間點製備）？
   → 有 → design formula 加 batch；PCA 第一主成份會被 batch 主導
```

### 資料品質紅線（啟動前用戶確認）

| 指標 | 建議確認方式 | 紅線 |
|------|------------|------|
| FASTQ 品質 | FastQC 或 fastp 快跑 | Q20 < 80% |
| 總 reads 數 | `zcat r1.fastq.gz \| wc -l` / 4 | < 5M reads/樣本 |
| 各樣本 reads 數差異 | 目測 FastQC summary | 最大/最小 > 10x |
| rRNA 污染率 | FastQC overrepresented sequences | > 30% rRNA |

### 進行中監控（STAR mapping 後）

```bash
# 跑完 STAR 後立即確認，不等 featureCounts
grep "Uniquely mapped reads %" runs/<run_id>/output/star_log/*_Log.final.out

# 紅線：任一樣本 < 50% → 停下來查，不要繼續
# 常見原因：genome/GTF 版本不匹配、物種錯誤、大量 adapter 殘留
```

---

## wgs_snp（WGS SNP/INDEL calling）

### 資料進入前必問

```
1. 樣本是什麼物種？二倍體還是多倍體？
   → ploidy 設定錯誤會讓所有 variant 頻率計算錯誤

2. 幾個樣本？預計 GPU 時間是否可接受？
   → 估算：樣本數 × 5h（30x WGS）= DeepVariant GPU 時間
   → 60 樣本 ≈ 12.5 天；啟動前讓用戶確認時程

3. 每個樣本的預期覆蓋深度（sequencing depth）？
   → < 10x：variant calling 信度低（DeepVariant 仍可跑，但 GQ 普遍低）
   → > 60x：確認 dv_shards 不要超過 CPU cores

4. 是否有參考基因組 + annotation（GTF/GFF）？
   → 無 annotation → SnpEff 步驟自動跳過（沒問題）
   → 有 annotation → 確認 GTF/GFF 版本與 genome 版本對應

5. 原始 FASTQ 是否已做過任何 QC？
   → 若已做過 adapter trimming → 告知 fastp 可以輕跑（--length_required 即可）
```

### 資料品質紅線

| 指標 | 紅線 | 後果 |
|------|------|------|
| FASTQ Q30 | < 70% | 大量 low-quality variant |
| 預估覆蓋深度 | < 5x/樣本 | min_depth filter 後幾乎沒有 PASS variant |
| 重複率（markdup 後） | > 40% | fixmate 排序錯誤（見 diagnose-snakemake-failure） |
| 重複率（markdup 後） | < 1% | markdup 可能未正確執行 |

> [!IMPORTANT]
> **Pilot test 必查重複率**：60 樣本全量跑完才發現 markdup 問題 = 損失 11 天。
> 2 樣本 × 最小 3 條 scaffold 的 pilot test 結束後，一定要看 `samtools flagstat`。

---

## genome_annotation（基因組注釋）

### 資料進入前必問（最複雜，問題最多）

```
1. 基因組大小是多少？
   → > 1.5 GB → 強烈建議 repeat_tool: "edta"（RepeatModeler 可能跑數週）
   → < 500 MB → repeat_tool: "repeatmodeler" 或 "both" 皆可

2. 是否有 RNA 證據？（最關鍵的問題）
   → 有 IsoSeq HiFi 資料？→ 設定 isoseq_dir（最佳品質）
   → 有短讀長 RNA-seq BAM？→ 設定 rnaseq_bam
   → 都沒有？→ 蛋白質僅版 BRAKER3，品質會明顯下降，需告知用戶

   教訓：RNA 資料已存在但未映射 → 等 3 小時映射 vs 重跑 BRAKER3 4 天
   永遠等 RNA 映射完成。

3. BUSCO lineage 是否選對？
   → 脊椎動物：vertebrata_odb10 或 actinopterygii_odb10
   → 軟體動物：mollusca_odb10
   → 昆蟲：insecta_odb10
   → BUSCO v5 vs v6：odb12 資料集只能搭配 BUSCO v6

4. 蛋白質資料庫（protein_fasta）是否為近緣物種？
   → 越近緣越好（同目/同科）
   → Metazoa.fa 是通用 fallback，但近緣物種效果更好

5. 是否有現成的 repeat library？
   → 有 → 可跳過 RepeatModeler/EDTA，直接指定 custom_repeat_library
   → 無 → 按步驟跑
```

### 資料品質紅線

| 指標 | 紅線 | 處置 |
|------|------|------|
| 基因組 N50 | < 100 kb | 組裝太碎，注釋品質低；建議先評估是否重新組裝 |
| BUSCO C%（組裝） | < 70% | 組裝本身有問題；注釋之前先確認原因 |
| FASTA header 數（rename 後） | 與原始不一致 | `\1` 轉義 bug（見 diagnose-snakemake-failure） |

---

## scrnaseq（單細胞 RNA-seq）

### 資料進入前必問

```
1. 10x Genomics chemistry 版本？
   → v2 / v3 / v3.1 → 影響 STARsolo 的 CB length 和 UMI length 設定

2. 預期細胞數？（實驗設計）
   → STARsolo cell calling 用 EmptyDrops_CR；結果應 ±20% 內
   → 差異 > 2x → 確認 chemistry 版本或 input cell loading 是否正確

3. 是否有 ambient RNA 污染的疑慮？
   → 通常 < 10% OK；> 20% 建議後處理（SoupX 或 DecontX）
   → Hermes 目前無自動 ambient correction → 需告知用戶

4. 多樣本是否需要整合（Harmony / Seurat integration）？
   → Hermes scrnaseq workflow 目前跑單一樣本 Leiden clustering
   → 多樣本整合目前需手動，或在 report workflow 後處理
```

### 資料品質紅線

| 指標 | 紅線 | 解讀 |
|------|------|------|
| Median genes/cell | < 200 | 死細胞比例高；考慮重新做 library |
| Median genes/cell | > 8000 | doublet 比例高（scDblFinder 會處理） |
| Mt% 中位數 | > 25% | 大量死細胞；結果不可信 |
| Valid barcode % | < 50% | chemistry 版本設定錯誤 |

---

## metagenome（HiFi 宏基因組）

### 資料進入前必問

```
1. HiFi reads 平均長度和 accuracy？
   → 建議 > 10 kb，accuracy > Q20（99%）
   → 短 reads（< 5 kb）→ hifiasm-meta 組裝效果差

2. 樣本來源（gut microbiome vs 環境）？
   → 宿主 DNA 比例？
   → 宿主 reads > 30% → 建議先做 host depletion（bowtie2 align to host genome）

3. 預估 MAG 數量？
   → 複雜環境（土壤、海水）→ MAG recovery 通常低（< 50 MAG）
   → 腸道 → 通常 100–500 MAG

4. 目標是 taxonomy 還是 MAG functional annotation？
   → taxonomy 只（Kraken2）→ 資料需求低，數小時可完成
   → MAG 重建（完整流程）→ 需要高覆蓋深度（> 50x 總 depth）
```

---

## pool_seq（Pool-seq Fst 掃描）

### 資料進入前必問

```
1. Pool 大小（每個 pool 幾個個體）？
   → Hudson Fst 計算需要 pool size 來校正；設定 pool_size 參數

2. 各 pool 間 coverage 是否均一？
   → 差異 > 5x → Fst 結果偏差；建議 downsampling 再分析

3. 是否有對應的 reference genome + annotation？
   → 無 annotation → 仍可跑 SNP + Fst scan，只是沒有基因層級注釋
```

---

## 快速決策：「這批資料可以啟動嗎？」

```
啟動前三問：
1. 資料完整性：gzip -t 沒有錯誤？
2. 樣本設計：replicate 數量足夠、分組合理？
3. 特定 workflow 的最關鍵問題（見上方各節）有回答？

三問都是 → 建議先跑 dry_run=True 預覽 DAG，用戶確認後啟動
任何一問有疑慮 → 先記錄疑慮，告知用戶風險，讓用戶明確確認後才啟動
```

> [!NOTE]
> 這個 checklist 是對話式把關，不是硬性阻擋。
> 如果用戶確認「我知道這個風險，繼續」→ 記錄在 run 的 config comment，然後啟動。
> 永遠讓用戶最終決定，小賀負責把風險說清楚。
