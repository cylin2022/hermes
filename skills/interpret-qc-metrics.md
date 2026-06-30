# interpret-qc-metrics

> **Skill**: 解讀各 workflow 的 QC 指標，判斷是否可進入下一步
> **觸發時機**: 用戶問「QC 結果如何？」「這個數字正常嗎？」「可以繼續跑嗎？」

---

## WGS / wgs_snp

### BAM 品質

| 指標 | 正常範圍 | 警示 | 失敗 |
|------|---------|------|------|
| 重複率（markdup） | 5–25% | 25–40% | > 40% 或 < 1% |
| Mapping rate | > 85% | 60–85% | < 60% |
| 平均覆蓋深度 | > 10x（分析）| 5–10x（低覆蓋） | < 5x |
| SO: tag（最終 BAM） | `coordinate` | — | `queryname`（排序錯誤） |

**重複率 > 40%**：`fixmate` 輸入未 name-sorted → 查 samtools 排序順序
**重複率 < 1%**：`markdup` 可能未正確執行 → 查 shell block

### VCF 品質

| 指標 | 正常範圍 | 備注 |
|------|---------|------|
| PASS variant 數（SNP） | 10k–5M（視基因組大小） | 非模式生物變異大，範圍廣 |
| Ti/Tv ratio | 2.0–2.5（WGS coding）| < 1.5 代表 false positive 高 |
| gVCF records/樣本 | 10M–100M | 小染色體測試可接受 10k–1M |
| GQ < 20 比例 | < 10% | > 30% 代表整體品質低 |

---

## RNA-seq / rnaseq

### STAR Mapping

| 指標 | 正常範圍 | 警示 | 失敗 |
|------|---------|------|------|
| Uniquely mapped % | > 70% | 50–70% | < 50% |
| Multi-mapped % | < 20% | 20–30% | > 30% |
| Unmapped: too short | < 10% | 10–20% | > 20% |

**Uniquely mapped < 50%**：
1. 先確認 reads 是否來自同一物種（污染？錯誤的參考基因組版本？）
2. 確認 genome 版本與 GTF 版本匹配
3. 考慮 `--outFilterMismatchNoverLmax 0.1`（對多型性高的非模式生物放寬）

### featureCounts 分配率

| 指標 | 正常範圍 | 警示 |
|------|---------|------|
| Assigned % | > 60% | < 40% |
| Unassigned_Nofeatures | < 30% | > 40% → strandedness 設定錯誤 |

**Unassigned_Nofeatures > 40%**：
- 執行 strandedness 檢查：`infer_experiment.py -r annotation.bed -i sample.bam`
- 結果 → featureCounts `-s` 參數對應：
  - 0.5/0.5 → `-s 0`（unstranded）
  - > 0.8 forward → `-s 1`
  - > 0.8 reverse → `-s 2`（最常見，Illumina TruSeq）

### DESeq2 結果

| 指標 | 正常 | 注意 |
|------|------|------|
| DE genes（padj < 0.05, |LFC| > 1） | 50–5000 | 0 → 統計power不足；> 10000 → 確認對照組正確 |
| Size factors | 0.5–2.0 | 極端值代表樣本間 library size 差異過大 |
| Dispersion estimate | 收斂 | 若發散 → 重複數不足或樣本品質差異大 |

---

## Genome Annotation / genome_annotation

### BUSCO

| 組別 | C% 正常範圍 | 解讀 |
|------|-----------|------|
| vertebrata_odb10 | > 85% | 脊椎動物參考 |
| metazoa_odb10 | > 70% | 廣用無脊椎動物 |
| 任意 lineage | S% > 80% | S = Single copy；D% > 5% 代表組裝問題 |

**C% < 60%**：
1. 先確認使用正確的 lineage（例如 `vertebrata` vs `actinopterygii`）
2. 若 lineage 正確 → 組裝品質本身有問題（N50 < 1 Mb？contig 太碎？）
3. 注意：BUSCO v5 不支援 odb12 資料集 → 改用 BUSCO v6

### RepeatModeler / EDTA

| 指標 | 正常範圍 | 備注 |
|------|---------|------|
| Repeat library 大小 | > 100 families | 小基因組 > 50 可接受 |
| Masked genome % | 20–80% | 魚類通常 30–50%；昆蟲 20–40% |

**Repeat library < 50 families**：
- 確認基因組大小（< 100 Mb 的小基因組正常）
- RepeatModeler 是否完成所有 Rounds（查 log 中 `Round` 數量）

### BRAKER3 Gene Model

| 指標 | 正常範圍 | 備注 |
|------|---------|------|
| 蛋白質編碼基因數 | 脊椎動物 15k–25k；軟體動物 20k–30k | |
| 平均外顯子數/基因 | 5–10 | |
| transcript_id 行數（braker.gtf） | > 50,000 | |

**基因數 < 5000**：
- 確認 BUSCO C% > 70%（先確認組裝品質）
- 確認有 RNA 證據（isoseq_dir 或 rnaseq_bam）
- 蛋白質僅版 BRAKER3 品質顯著較差 → 補齊 RNA 再重跑

---

## scRNA-seq / scrnaseq

### STARsolo Mapping

| 指標 | 正常範圍 | 備注 |
|------|---------|------|
| Reads with valid CB+UMI | > 60% | CB = Cell Barcode |
| Estimated number of cells | 接近實驗預期值 ±20% | |

### Scanpy / Leiden Clustering

| 指標 | 正常 | 注意 |
|------|------|------|
| Median genes/cell | 500–5000 | < 200 → 死細胞比例高；> 8000 → doublet |
| Mitochondria % | < 20% | > 30% → 大量死細胞 |
| n_obs（細胞數） | 接近預期 | |

**doublet 比例 > 10%（scDblFinder）**：
- 確認上游 CB filtering 閾值
- 建議重新以更嚴格的 STARsolo CB min_reads 過濾

---

## 快速判斷：「可以繼續嗎？」

```
一律先確認這三項：
1. Snakemake 最後一行是否 "N of N steps (100%) done"？
2. 所有輸出檔案 owner 是 cylin（不是 root）？
3. 品質指標是否在上表的「正常範圍」？

三項都是 → 可繼續
任何一項否 → 先診斷再繼續（參考 diagnose-snakemake-failure skill）
```

---

## MultiQC 報告位置

每個 workflow 完成後，MultiQC 報告在：
```
runs/<run_id>/output/multiqc/multiqc_report.html
```

用 `read_file(path)` 讀取 summary 表格，或讓用戶直接在瀏覽器開啟。
