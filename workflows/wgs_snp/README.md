# WGS SNP/INDEL Calling — Non-Model Organisms

Whole-genome SNP and INDEL discovery from short-read WGS data, using **DeepVariant + GLnexus**,
optimised for non-model organisms.

## Why DeepVariant over GATK HaplotypeCaller?

| | GATK HaplotypeCaller | DeepVariant |
|--|---|---|
| Algorithm | Heuristic rules (human-tuned) | Neural network (generalises to new genomes) |
| BQSR required | Yes (impossible without known-variant DB) | **No** |
| Filter tuning | Hard-filter thresholds need manual adjustment | **No — GQ ≥ 20, DP ≥ 5 is sufficient** |
| SNP accuracy | High | **Higher** (better Precision/Recall, Ts/Tv) |
| INDEL accuracy | Good | Comparable, improving |
| Joint calling | Native GVCF → GenotypeGVCFs | gVCF → **GLnexus** (10× faster) |

## Pipeline steps

1. **fastp** — adapter trimming, quality filtering, per-sample HTML QC report
2. **MultiQC** — aggregate QC across all samples
3. **BWA-MEM2** — fast alignment with read-group tags
4. **samtools fixmate + sort + markdup** — coordinate sort and PCR duplicate marking
5. **DeepVariant** — per-sample variant calling in gVCF mode (`--model_type WGS`)
   - Runs via **Docker `--gpus all`** on NVIDIA RTX Ada 6000 (CUDA 12); ~30 min/sample
   - Falls back to conda CPU package when `use_gpu: false`
6. **GLnexus** — joint genotyping across all samples (`--config DeepVariant`)
7. **bcftools filter** — apply GQ ≥ 20 and DP ≥ 5 soft-filters, then extract PASS SNPs/INDELs
8. **SNPeff** — functional annotation with auto-built custom database from genome + GTF/GFF
9. **bcftools stats** — Ts/Tv ratio, per-sample SNP counts, summary statistics

## Key outputs

| File | Description |
|------|-------------|
| `vcf/filtered/snps.PASS.vcf.gz` | PASS SNPs (GQ ≥ min_gq, DP ≥ min_depth) |
| `vcf/filtered/indels.PASS.vcf.gz` | PASS INDELs |
| `annotation/snps.annotated.vcf.gz` | SNPs with functional ANN field from SNPeff |
| `annotation/snpeff_summary.html` | SNPeff annotation summary report |
| `stats/variant_summary.txt` | SNP/INDEL counts, Ts/Tv, per-sample stats |
| `alignment/{sample}.markdup.bam` | Final alignments (duplicate-marked) |
| `multiqc/multiqc_report.html` | QC report |

## Quality filtering philosophy

DeepVariant was trained to output well-calibrated probabilities. The FILTER field already
contains `PASS` vs `RefCall` for high-quality calls. The additional `GQ` and `DP` filters
catch edge cases (very low-coverage sites, uncertain genotypes). Unlike GATK, you do **not**
need to tune QD, FS, MQ, ReadPosRankSum thresholds.

## SNPeff custom database

The pipeline automatically builds a SNPeff database:
- `species_name` → database name (alphanumeric + underscores only)
- GTF preferred; GFF3 supported as fallback
- Set both `gtf` and `gff` to `""` to skip annotation

## Notes on ploidy

- `ploidy: 2` → diploid (stable, recommended)
- `ploidy: 4+` → experimental support in DeepVariant v1.6; test with dry_run first
- For strict polyploid variant calling, consider GATK HC as an alternative

## Estimated runtime (160 cores, 10 samples, 1 Gb genome, 30× coverage)

| Step | Time (10 samples, GPU) | Time (10 samples, CPU) |
|------|------------------------|------------------------|
| fastp | ~30 min | ~30 min |
| BWA-MEM2 + markdup (per sample) | ~1–2 h | ~1–2 h |
| DeepVariant (per sample) | **~30 min** | ~3–4 h |
| GLnexus joint calling | ~15 min | ~15 min |
| Filtering + annotation | ~20 min | ~20 min |
| **Total** | **~2–4 h** | **~6–12 h** |

## GPU setup (RTX Ada 6000 / CUDA 12)

The workflow uses Docker + `--gpus all` for DeepVariant (default: `use_gpu: true`).
Prerequisites already confirmed on this machine:
- Docker 29.1.3 ✓
- CUDA 12.0 ✓
- NVIDIA RTX 6000 Ada 49 GB ✓

`nvidia-container-toolkit` must be installed for `--gpus all` to work:
```bash
# Check
docker run --rm --gpus all nvidia/cuda:12.0-base-ubuntu22.04 nvidia-smi
# Install if missing
sudo apt install nvidia-container-toolkit && sudo systemctl restart docker
```

To use CPU instead: set `use_gpu: false` in config.
