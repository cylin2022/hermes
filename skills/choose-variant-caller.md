# choose-variant-caller

> **Skill**: 根據實驗設計選擇正確的 variant caller 及流程配置
> **觸發時機**: 用戶描述 variant calling 需求時

---

## 決策樹：選擇 Variant Caller

```
問題一：樣本是否為模式生物（人、鼠、斑馬魚）？
│
├── 是 → 考慮 GATK HaplotypeCaller（有成熟 VQSR filter model）
│         但 Hermes 目前未整合 GATK，建議仍使用 DeepVariant + GLnexus
│
└── 否（非模式生物，如吳郭魚、文蛤、蝦等）→ 繼續下方

問題二：資料類型？
│
├── WGS（全基因組，多個個體）
│   └── → wgs_snp workflow（DeepVariant + GLnexus）
│
├── Pool-seq（混合樣本，無個體資訊）
│   └── → pool_seq workflow（BWA-MEM2 + bcftools + Hudson Fst）
│
├── RNA-seq（轉錄組 variant）
│   └── → 目前無專用 workflow；需手動或 rnaseq workflow 後處理
│
└── Amplicon / Target sequencing
    └── → 目前無專用 workflow；建議諮詢需求
```

---

## 各 Caller 特性速查

| 工具 | 適用場景 | 本機狀態 | 備注 |
|------|---------|---------|------|
| **DeepVariant 1.10.0-gpu** | WGS、非模式生物 | ✅ 已整合（wgs_snp） | GPU 加速，準確率業界最高 |
| **GLnexus** | 多樣本 joint calling | ✅ 已整合（wgs_snp） | 接收 DeepVariant gVCF |
| **bcftools call** | Pool-seq、快速 SNP | ✅ 已整合（pool_seq） | 不適合個體 gVCF |
| **GATK HaplotypeCaller** | 模式生物 | ❌ 未整合 | 需要 VQSR 訓練集 |
| **Freebayes** | 小型非模式生物計畫 | ❌ 未整合 | 無 GPU 加速 |

---

## wgs_snp Workflow 配置要點

### 何時用 GPU（use_gpu: true）

- 預設應**開啟 GPU**（`use_gpu: true`，`dv_docker_image: google/deepvariant:1.10.0-gpu`）
- GPU 加速 DeepVariant 約快 10–20 倍
- GPU 資源設定：Snakemake 啟動時加 `--resources gpu=1`

```yaml
# config.yaml 範例
use_gpu: true
dv_docker_image: "google/deepvariant:1.10.0-gpu"
dv_shards: 64       # 建議 = CPU threads；GPU 機器可設高
threads: 64
mem_gb: 128
```

### 樣本數與資源評估

| 樣本數 | 建議 dv_shards | 預估 DeepVariant 時間（/樣本） | 總 GPU 時間 |
|--------|--------------|-------------------------------|------------|
| 1–10   | 32–64        | 3–6 h（WGS ~30x）            | < 2 天     |
| 11–60  | 64           | 3–6 h（WGS ~30x）            | 7–15 天    |
| > 60   | 64           | 建議分批，或確認 GPU 供應充足 | > 15 天    |

> [!IMPORTANT]
> DeepVariant 為序列執行（1 GPU 一次 1 個樣本）。
> 60 個樣本的 GPU 時間估算：60 × 5h = 300h ≈ 12.5 天。
> 啟動前請確認這個時程可接受。

### Ploidy 設定

```yaml
ploidy: 2    # 二倍體（大多數動植物）
ploidy: 1    # 單倍體（雄性蜂、某些藻類）
ploidy: 4    # 四倍體（小麥、草莓等多倍體）
```

> 注意：GLnexus joint calling 的 ploidy 需與 DeepVariant 設定一致。

---

## Pool-seq vs WGS：選擇標準

| 情境 | 推薦 |
|------|------|
| 有個體 DNA，要找每個個體的 genotype | **wgs_snp**（DeepVariant + GLnexus） |
| 混池 DNA，只要群體頻率 + Fst | **pool_seq**（bcftools + Hudson Fst） |
| 有個體 DNA，但樣本數 > 200，資源有限 | 考慮分批 wgs_snp 或改 pool_seq |
| 要做 GWAS | **wgs_snp** → **snp_association** |
| 要做 genomic prediction | **wgs_snp** → **genomic_prediction** |

---

## SNP Filter 參數建議（非模式生物）

```yaml
# wgs_snp config_template.yaml
min_gq: 20         # Genotype Quality；< 20 視為低品質
min_depth: 5       # 最小覆蓋深度；< 5 信度不足
```

**下游 GWAS 前的額外 filter（在 snp_association workflow 中設定）**：
- MAF > 0.05（GEMMA LMM）
- Missing rate < 0.1（PLINK2 `--geno 0.1`）
- HWE p > 1e-6（`--hwe 1e-6`）

---

## SNP Annotation 配置（SnpEff）

wgs_snp workflow 使用 SnpEff 建立自定義資料庫。

**輸入需求**：
- `genome_fasta`：參考基因組
- `gtf` 或 `gff`：基因 annotation

**常見問題**：
- 只接受 `genes.gff`，不接受 `genes.gff3` → 副檔名必須用 `gff`
- NCBI GTF 缺少 CDS sequence → 加 `-noCheckCds -noCheckProtein`
- SnpEff 版本 > 5.2 需要 Java 17

**無 annotation 的情況**：
- `gtf` 和 `gff` 都留空 → SnpEff annotation rule 自動跳過
- 輸出 VCF 可用於後續 GWAS，只是沒有功能注釋

---

## 流程串接建議

```
wgs_snp
  └─ snp_association（GWAS：PLINK2 + GEMMA LMM + Fst scan）
       └─ 需要：filtered VCF + phenotype CSV

wgs_snp
  └─ genomic_prediction（ML 預測：GBLUP + LASSO + RF + XGBoost）
       └─ 需要：filtered VCF + phenotype CSV

wgs_snp + rnaseq
  └─ 整合 eQTL 分析（目前無 workflow，需手動）
```
