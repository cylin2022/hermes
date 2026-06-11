# SNP Association Analysis — Binary Phenotype + Related Individuals

Identifies SNPs associated with salt tolerance in Taiwan tilapia (30 tolerant vs 30 intolerant).
Uses two complementary strategies that reinforce each other:

| Strategy | Tool | Strength |
|----------|------|----------|
| **Mixed-model GWAS** | GEMMA LMM | Controls kinship; per-SNP significance test |
| **Fst outlier scan** | VCFtools | Detects genomic regions under divergent selection |

SNPs significant in both strategies are the highest-confidence candidates.

## Handling body-size confounding

The salt-tolerant group had larger body size at experiment start. Without correction,
any GWAS would conflate growth QTLs with salt-tolerance QTLs. This pipeline runs
**two GWAS in parallel** and classifies every SNP:

| Category | Definition | Interpretation |
|----------|-----------|----------------|
| `salt_specific` | p_salt < threshold, p_size ≥ threshold | True salt-tolerance SNP |
| `salt_fst` | salt_specific + Fst outlier | Strongest candidates |
| `pleiotropic` | Both p_salt and p_size significant | Same gene affects growth AND salt tolerance |
| `size_specific` | p_size < threshold, p_salt ≥ threshold | Growth QTL (confounder, filtered out) |
| `fst_only` | Fst outlier only | Population divergence without individual-level association |

**Key design**: GEMMA GWAS A includes body weight + length as continuous covariates (`-c`),
mathematically removing their effects before testing salt-tolerance association.

## Pipeline steps

1. **Prepare inputs** — parse phenotype + covariate files (body weight/length z-score normalised)
2. **VCF → PLINK QC** — MAF ≥ 5%, missingness, HWE; `--allow-extra-chr` for tilapia scaffolds
3. **LD pruning** — r² < 0.1 in 50-SNP windows for PCA/kinship
4. **PCA** (PLINK2) — visualise population structure
5. **Kinship matrix** (GEMMA -gk 1) — 60×60 centred GRM; absorbs family-structure inflation
6. **GWAS A** (GEMMA -lmm 4 + `-c` covariates) — salt tolerance, adjusted for body size
7. **GWAS B** (GEMMA -lmm 4) — body weight as phenotype; identifies growth QTLs
8. **Fst per-SNP + sliding window** (VCFtools)
9. **Plots** (R/CMplot) — Manhattan (adjusted GWAS + Fst), QQ + λ, PCA
10. **Annotate + classify** — three-way SNP classification + nearest gene (BEDTools)

## Key outputs

| File | Description |
|------|-------------|
| `plots/manhattan_gwas.pdf` | GWAS Manhattan plot (red line: p < 5×10⁻⁸; blue dashed: p < 1×10⁻⁵) |
| `plots/manhattan_fst.pdf` | Genome-wide Fst scan |
| `plots/qq_plot.pdf` | QQ plot + genomic inflation factor λ |
| `plots/pca.pdf` | PCA coloured by phenotype group |
| `annotation/candidate_snps.csv` | Candidates with p_lrt, Fst, beta, nearest gene |
| `annotation/candidate_genes.csv` | Unique candidate genes ranked by best p-value |
| `gwas/gwas.assoc.txt` | Full GEMMA output for all SNPs |
| `fst/persnp.weir.fst` | Per-SNP Fst for all SNPs |
| `fst/window.windowed.weir.fst` | Sliding-window Fst |

## Important notes for 10× WGS coverage

This workflow was designed with **10× per-individual coverage** in mind. Implications:

- **DeepVariant** in the `wgs_snp` workflow handles 10× well; use `min_depth: 3` there
- **Genotype error rate** is higher at 10× (~1–2%) vs 30× — use `maf ≥ 0.05` to reduce false associations from rare-allele miscalls
- **For improved accuracy**: consider running `ANGSD` with genotype likelihoods instead of
  hard-called genotypes; especially beneficial for PCA (PCAngsd) and population genetics stats
- **Statistical power**: with n=60 at 10×, detectable effect sizes are large (OR > 3).
  The extreme phenotype design (top vs bottom individuals) partially compensates for small n.

## Interpreting λ (genomic inflation factor)

| λ | Interpretation |
|---|---------------|
| 1.00–1.05 | Well-controlled; no inflation |
| 1.05–1.10 | Mild inflation; check PCA for stratification |
| > 1.10 | Significant inflation; consider adding PCs as covariates in GEMMA |

GEMMA's kinship matrix should absorb most family-structure inflation. If λ > 1.10,
increase kinship precision by using more LD-pruned SNPs.

## Workflow order

```
wgs_snp (60 samples) → snp.PASS.vcf.gz
          ↓
snp_association → candidate_snps.csv + candidate_genes.csv
```

## Estimated runtime (160 cores, ~500k QC-passing SNPs, 60 samples)

| Step | Time |
|------|------|
| VCF → PLINK QC | ~10 min |
| Kinship + PCA | ~20 min |
| GWAS (GEMMA) | ~1–2 h |
| Fst (per-SNP + window) | ~30 min |
| Plots + annotation | ~20 min |
| **Total** | **~2–3 h** |
