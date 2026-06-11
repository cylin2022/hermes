# Hermes Workflow Performance & Efficiency Audit

**Date:** 2026-06-12  
**Hardware:** 160 cores · 2.2 TiB RAM · NVIDIA RTX 6000 Ada (49 GB VRAM) · Ubuntu 24.04 · 66 TB fast storage  
**Auditor:** Claude Code (claude-sonnet-4-6)

---

## 1. Executive Summary

| Workflow | Snakefile Present | Overall Status | Primary Issues |
|---|---|---|---|
| `wgs_snp` | Yes | **Needs Tuning** | threads=16 default (10× underuse), GPU DV runs sequentially, markdup BAM not `temp()`, GLnexus missing `--mem-gbytes`, snpEff missing `-Xmx` |
| `snp_association` | Yes | **Needs Tuning** | No `threads:` on any rule, GEMMA single-threaded by design but PLINK2/VCFtools unused parallelism, undefined `PHENOTYPE` variable (runtime crash), `gwas_bodysize` always runs |
| `scrnaseq` | Yes | **Needs Tuning** | threads=16 default, STARsolo BAM not `temp()`, scDblFinder uses `SerialParam()` (single-core), `mem_gb=128` tight for STAR index build on large genomes |
| `rnaseq` | **No** (Snakefile missing) | **Incomplete** | Only config_template.yaml + README exist; no Snakefile |
| `atacseq` | **No** (Snakefile missing) | **Incomplete** | Only samplesheet_template.csv exists; no Snakefile or config |
| `genome_annotation` | **No** (Snakefile missing) | **Incomplete** | Only samplesheet_template.csv exists; no Snakefile or config |

**Score legend:** Ready = production-ready as-is · Needs Tuning = functional but leaves significant performance on the table · Incomplete = Snakefile not written

---

## 2. Per-Workflow Findings

---

### 2.1 `wgs_snp`

#### A. Thread/Core Allocation

**Finding:** The global `THREADS` variable defaults to `16` in config and is used verbatim for `bwa_mem2`, `deepvariant`, and `glnexus`. On 160 cores the server can run many samples in parallel, but each individual BWA-MEM2 or DeepVariant job only uses 16 of those cores even if it is the only job running.

**Config default (line 33):**
```yaml
threads: 16
```

**Recommended change** — raise the default to 64 so a single-sample run saturates a quarter of the machine, while still allowing 2-3 concurrent samples:
```yaml
threads: 64
```

**Finding:** `rule fastp` hardcodes `threads: 8` (line 76). fastp scales well to 16+ threads for large genomes.

```python
# BEFORE (line 76)
threads: 8

# AFTER
threads: 16
```

**Finding:** `rule sort_markdup` hardcodes `threads: 4` (line 141) and `mem_mb = 16 * 1024`. For 160-core hardware this is a bottleneck — samtools sort benefits substantially from more threads when sorting large BAMs.

```python
# BEFORE (lines 141-143)
threads: 4
resources:
    mem_mb = 16 * 1024,

# AFTER
threads: 8
resources:
    mem_mb = 32 * 1024,
```

**Finding:** `rule flagstat` has no `threads:` block; it is single-threaded. `samtools flagstat` accepts `-@` for I/O threads. Not critical but free wins:

```python
# ADD to rule flagstat
threads: 4
shell: "samtools flagstat -@ {threads} {input.bam} > {output}"
```

#### B. Memory Allocation

**Finding:** Default `mem_gb: 64` in config. For DeepVariant with `--num_shards=64`, memory peaks at ~3–4 GB per shard, so 64 threads × 4 GB = ~256 GB needed at peak. At `mem_gb=64` the container will OOM if `threads` is raised to 64 without also raising `mem_gb`.

**Recommended config change:**
```yaml
threads: 64
mem_gb: 256   # 64 shards × ~4 GB/shard peak
```

**Finding:** `rule glnexus` uses `mem_mb = MEM_GB * 1024` but passes no explicit memory cap to `glnexus_cli`. GLnexus defaults to using all available RAM, which is correct, but the shell command should add `--mem-gbytes` so the resource request is enforced:

```python
# BEFORE (lines 263-273)
shell:
    """
    rm -rf {params.db}
    mkdir -p $(dirname {output.vcf})
    glnexus_cli \
        --config DeepVariant \
        --dir {params.db} \
        --threads {threads} \
        {input.gvcfs} | \
    bcftools view --output-type z --output {output.vcf}
    tabix -p vcf {output.vcf}
    """

# AFTER — add mem cap and bcftools thread flag
shell:
    """
    rm -rf {params.db}
    mkdir -p $(dirname {output.vcf})
    glnexus_cli \
        --config DeepVariant \
        --dir {params.db} \
        --threads {threads} \
        --mem-gbytes {MEM_GB} \
        {input.gvcfs} | \
    bcftools view --output-type z --threads 4 --output {output.vcf}
    tabix -p vcf {output.vcf}
    """
```

**Finding:** `rule snpeff_annotate` runs `snpEff ann` with no JVM heap flag. For large genomes snpEff will use default 256 MB heap and GC-thrash. Add `-Xmx`:

```python
# BEFORE (line 364)
snpEff ann \

# AFTER
snpEff ann -Xmx16g \
```

#### C. GPU Utilization

**Finding:** GPU is used correctly via `docker run --rm --gpus all`. The DeepVariant Docker image is `google/deepvariant:1.6.1` which includes CUDA support for the Ada RTX 6000. This is well-configured.

**Finding:** The GPU-enabled DeepVariant rule runs one sample at a time (sequential expand in `rule all`). With a 49 GB GPU and DeepVariant using ~8–12 GB VRAM per sample, it is theoretically possible to run 2 concurrent samples. However, Snakemake's resource system does not natively model VRAM, so this would require a `resources: gpu=1` constraint to prevent two simultaneous GPU jobs. The current design is safe but leaves ~50% GPU idle during each DV run.

**Recommended addition** — add GPU resource tracking to DeepVariant rule and pass `--resources gpu=1` at launch:

```python
# ADD to both GPU deepvariant rule resource blocks
resources:
    mem_mb = MEM_GB * 1024,
    gpu    = 1,               # sentinel: prevents >1 concurrent DV GPU job
```

Then launch Snakemake with `--resources gpu=1`.

#### D. I/O and Temp Files

**Critical finding:** `rule bwa_mem2` correctly marks the namesorted BAM as `temp()` (line 122). However, `rule sort_markdup` output `{sample}.markdup.bam` and `{sample}.markdup.bam.bai` (lines 139-140) are **not** marked `temp()`. After DeepVariant produces the gVCF, these large per-sample BAMs (~20–100 GB each) are never cleaned up. For 30+ samples this is hundreds of GB of stranded data.

```python
# BEFORE (lines 139-140)
output:
    bam = str(OUTDIR / "alignment" / "{sample}.markdup.bam"),
    bai = str(OUTDIR / "alignment" / "{sample}.markdup.bam.bai"),

# AFTER — mark BAM as temp but keep BAI for flagstat; better: keep both and clean after DV
# Option A: mark both temp (user loses BAMs after pipeline)
output:
    bam = temp(str(OUTDIR / "alignment" / "{sample}.markdup.bam")),
    bai = temp(str(OUTDIR / "alignment" / "{sample}.markdup.bam.bai")),

# Option B (recommended): add a config flag keep_bam: false and conditionally wrap
```

**Finding:** DeepVariant's `{sample}_tmp` intermediate directory is already cleaned (`rm -rf {params.tmp_dir}`) inside the shell block. Correct.

**Finding:** `rule filter_variants` marks `all.filtered.vcf.gz` as `temp()` (line 282). Correct.

#### E. Tool Version / Best Practice Flags

**Finding:** `bcftools view` in `rule glnexus` (line 271) does not use `--threads`. Add `--threads 4` for faster VCF compression.

**Finding:** `snpEff build` and `snpEff ann` do not pass `-Xmx` (JVM heap). For genomes >500 MB, snpEff will be slow or crash with default 256 MB.

**Finding:** `bwa-mem2` version pinned to 2.2.1 (released 2021). Version 2.2.1 is stable and current; no change needed.

**Finding:** `samtools fixmate` in the BWA pipeline uses `-@ 4` hardcoded (line 132). Should use `{threads}` or at minimum a higher value since fixmate is I/O-bound:

```python
# BEFORE (line 132)
samtools fixmate -m -@ 4 - {output.bam}

# AFTER
samtools fixmate -m -@ {threads} - {output.bam}
```

#### F. Missing Features / Gaps

**Finding:** No VQSR or equivalent post-filter confidence score. DeepVariant's calibrated QUAL scores make VQSR unnecessary, so the current GQ+DP filter is appropriate. No gap here.

**Finding:** No contamination check (e.g., VerifyBamID) or ancestry outlier detection step before variant calling. Not critical but recommended for production.

**Finding:** No `--rerun-incomplete` flag documentation or Snakemake profile. Add to README or MCP runner.

**CRITICAL BUG:** The `snpeff_annotate` rule only annotates SNPs (`snps.PASS.vcf.gz`) but not INDELs. The INDEL annotation output is declared in `rule all` via `if ANNOT_FILE` but only `snps.annotated.vcf.gz` is produced. Add a corresponding `rule snpeff_annotate_indels` or annotate the merged filtered VCF instead.

---

### 2.2 `snp_association`

#### A. Thread/Core Allocation

**Critical finding:** No single rule in this Snakefile has a `threads:` directive. All rules run single-threaded. PLINK2 and VCFtools both support parallelism.

```python
# ADD threads to rule vcf_to_plink:
threads: 16
shell:
    """
    plink2 \
        --vcf {input.vcf} \
        --threads {threads} \
        ...
    """

# ADD threads to rule ld_prune:
threads: 16
shell:
    "plink2 --bfile {params.p} --threads {threads} --indep-pairwise 50 10 0.1 ..."

# ADD threads to rule pca:
threads: 16
shell:
    "plink2 --bfile {params.p} --threads {threads} --extract {input.keep} --pca {params.n} ..."
```

**Finding:** GEMMA (`rule kinship`, `rule gwas_salt_adjusted`, `rule gwas_bodysize`) is inherently single-threaded in its LMM implementation (gemma 0.98.5 does not parallelize the BFGS solver). No change possible there.

**Finding:** VCFtools (`rule fst_persnp`, `rule fst_window`) is single-threaded by design. The computation is fast (~minutes for 60 samples). No change needed.

#### B. Memory Allocation

**Finding:** No `resources: mem_mb` blocks on any rule. PLINK2's VCF-to-BED conversion for a large cohort can use 8–16 GB. GEMMA kinship computation scales as O(N² × M) where M = SNPs; for 60 samples × 5M SNPs this peaks at ~20–30 GB.

**Add resources blocks:**
```python
rule vcf_to_plink:
    resources:
        mem_mb = 16 * 1024,

rule kinship:
    resources:
        mem_mb = 32 * 1024,

rule gwas_salt_adjusted:
    resources:
        mem_mb = 32 * 1024,

rule gwas_bodysize:
    resources:
        mem_mb = 32 * 1024,
```

#### C. GPU Utilization

No GPU-acceleratable steps in this workflow. GEMMA does not support GPU acceleration. No gap.

#### D. I/O and Temp Files

**Finding:** PLINK binary files (`.bed`, `.bim`, `.fam`) and the LD-pruned subset files are large intermediates that are never cleaned. After the GWAS is complete they are not needed. Mark as `temp()`:

```python
rule vcf_to_plink:
    output:
        bed = temp(PFX + ".bed"),
        bim = temp(PFX + ".bim"),
        fam = temp(PFX + ".fam"),

rule ld_prune:
    output:
        keep = temp(str(OUTDIR / "plink" / "pruned.prune.in")),
        excl = temp(str(OUTDIR / "plink" / "pruned.prune.out")),
```

Note: if users want to keep PLINK binaries for downstream tools, add a `keep_plink: true` config switch.

#### E. Tool Version / Best Practice Flags

**CRITICAL BUG — NameError at parse time:** Line 281 references `PHENOTYPE` which is never defined in the Snakefile:

```python
# Line 281 (rule plots, input block):
pheno = PHENOTYPE,   # <-- NameError: name 'PHENOTYPE' is not defined
```

The variable should be `METADATA` (the metadata CSV path):

```python
# BEFORE
pheno = PHENOTYPE,

# AFTER
pheno = METADATA,
```

This bug will cause Snakemake to crash at DAG-building time, before any rules run.

**Finding:** `rule gwas_bodysize` is unconditionally included in the DAG even when `HAS_COV = False`. The `rule all` input correctly gates on `HAS_COV`, but the rule body still tries to use `input.pheno` which will be an empty file written by `prepare_inputs`. GEMMA will produce an empty or malformed output. The rule should be conditionally defined:

```python
# AFTER rule gwas_salt_adjusted, wrap gwas_bodysize:
if HAS_COV:
    rule gwas_bodysize:
        ...
```

**Finding:** `plink2 --hwe {params.hwe} keep-controls` — the `keep-controls` flag is space-separated from `--hwe`, which is a PLINK2 modifier syntax. This is correct PLINK2 syntax; no change needed.

#### F. Missing Features / Gaps

**Finding:** No `resources:` mem blocks anywhere (addressed above).

**Finding:** No population stratification check output beyond PCA eigenvectors (no scree plot, no population outlier flagging). Not critical for n=60 but worth adding.

**Finding:** The `manhattan_plot.R` script reads input via `snakemake@input[["gwas"]]` (line 9) but the `rule plots` Snakemake input is named `gwas_salt` (line 276), not `gwas`. This mismatch will cause an R-level error when the script tries to read the file.

```r
# BEFORE (manhattan_plot.R line 9):
gwas_file   <- snakemake@input[["gwas"]]

# AFTER:
gwas_file   <- snakemake@input[["gwas_salt"]]
```

---

### 2.3 `scrnaseq`

#### A. Thread/Core Allocation

**Finding:** Default `threads: 16` in config. Same structural issue as wgs_snp — should be raised for single-run use.

```yaml
# config_template.yaml recommended change:
threads: 32   # good balance for STARsolo; Scanpy PCA benefits from more
```

**Finding:** `rule fastqc` hardcodes `threads: 4`. FastQC parallelism scales by file pair, not intra-file. For 2 files this is already optimal. No change needed.

**Finding:** `rule doublet_detection` (scDblFinder) uses `SerialParam()` in the R script (line 18 of doublet_detection.R). This explicitly disables parallelism. For samples with 10,000+ cells, scDblFinder can take several minutes per sample single-threaded. Switch to `MulticoreParam`:

```r
# BEFORE (doublet_detection.R line 18):
sce <- scDblFinder(sce, BPPARAM = SerialParam())

# AFTER:
nworkers <- if (!is.null(snakemake@threads)) snakemake@threads else 4L
sce <- scDblFinder(sce, BPPARAM = MulticoreParam(workers = nworkers))
```

And add `threads: 8` to `rule doublet_detection` in the Snakefile:

```python
rule doublet_detection:
    ...
    threads: 8
    conda: "envs/r_scrna.yaml"
    script: "scripts/doublet_detection.R"
```

**Finding:** `rule scanpy_pipeline` correctly passes `snakemake.threads` as `n_jobs` and sets `sc.settings.n_jobs`. However, Scanpy's `rank_genes_groups` with `method="wilcoxon"` does not parallelize via `n_jobs`. The PCA (`svd_solver="arpack"`) is also single-threaded in scipy. For this hardware, consider using `svd_solver="randomized"` for speed at slight accuracy tradeoff, or RAPIDS GPU-accelerated Scanpy (`rapids-singlecell`) for the PCA/UMAP steps.

#### B. Memory Allocation

**Finding:** `mem_gb: 128` in config. STAR genome index build for a mammalian genome (~3 GB genome) needs ~33 GB RAM. The subsequent STARsolo alignment against a loaded index uses 33 GB + read buffer. 128 GB is sufficient but conservative. For very large genomes (e.g., wheat at 17 GB) this will OOM; the config comment should warn users.

**Finding:** Scanpy `rank_genes_groups` with Wilcoxon on 50,000+ cells × 20,000 genes can use 30–60 GB. `mem_gb=128` provides headroom but should be documented.

**Finding:** No `resources: mem_mb` block on `rule doublet_detection`. scDblFinder for 30,000 cells uses ~8–16 GB.

```python
rule doublet_detection:
    resources:
        mem_mb = 16 * 1024,
```

#### C. GPU Utilization

**Finding:** Scanpy uses CPU for PCA, neighbor graph, and UMAP. The RTX 6000 Ada (49 GB VRAM) could dramatically accelerate all three steps via `rapids-singlecell` (cuML/cuGraph backend). For 50,000+ cells, GPU UMAP is 10–50× faster. This is an optional but high-value enhancement.

**Recommendation:** Add an optional `use_gpu: false` config flag and a `envs/rapids_scanpy.yaml` environment:

```yaml
# envs/rapids_scanpy.yaml (new file)
channels:
  - rapidsai
  - nvidia
  - conda-forge
dependencies:
  - rapids-singlecell>=0.10
  - cuda-version=12.*
```

And in `scanpy_pipeline.py`:

```python
# Near top of script, after imports:
use_gpu = snakemake.config.get("use_gpu", False)
if use_gpu:
    import rapids_singlecell as rsc
    rsc.get.anndata_to_GPU(adata)
    # Replace sc.tl.pca, sc.pp.neighbors, sc.tl.umap with rsc equivalents
```

#### D. I/O and Temp Files

**Critical finding:** `rule starsolo` produces `{sample}/Aligned.sortedByCoord.out.bam` (line 127). For a typical scRNA-seq sample this BAM is 5–30 GB. After `rule doublet_detection` completes, this BAM is never used again in the pipeline. It should be marked `temp()`:

```python
# BEFORE (line 127):
bam = str(OUTDIR / "starsolo" / "{sample}" / "Aligned.sortedByCoord.out.bam"),

# AFTER:
bam = temp(str(OUTDIR / "starsolo" / "{sample}" / "Aligned.sortedByCoord.out.bam")),
```

For 10 samples × 15 GB = 150 GB saved automatically.

**Finding:** STARsolo genome-loaded shared memory (`--genomeLoad LoadAndRemove`) is not used. For running multiple samples sequentially, loading the genome once saves ~5–10 minutes per additional sample. Add `--genomeLoad LoadAndRemove` to `rule starsolo` shell (only when samples run sequentially, which is the Snakemake default with shared STAR index dependency).

#### E. Tool Version / Best Practice Flags

**Finding:** STAR 2.7.11a is current and appropriate.

**Finding:** `rule starsolo` does not pass `--outBAMsortingBinsN` to STAR. For large samples STAR's default bin count (50) is too low for sorting, causing temp file explosion. Add:

```python
# ADD to STAR --runMode alignReads call:
--outBAMsortingBinsN 100 \
--outBAMsortingThreadN {threads} \
```

**Finding:** scDblFinder R package version not pinned in `envs/r_scrna.yaml`. Add version pin for reproducibility:

```yaml
# envs/r_scrna.yaml
- bioconductor-scdblfinder=1.16
```

#### F. Missing Features / Gaps

**Finding:** No ambient RNA correction step (SoupX or CellBender) before doublet detection. For low-quality samples with high ambient contamination this can inflate false markers. Not critical but a best-practice gap.

**Finding:** No batch correction step (Harmony, scVI) in `scanpy_pipeline.py` even though `batch_key` is passed to HVG selection. The PCA and neighbors graph do not use any integration method. For multi-sample datasets with batch effects, Harmony should be added after PCA:

```python
# ADD after sc.tl.pca (scanpy_pipeline.py):
if len(samples) > 1:
    sc.external.pp.harmony_integrate(adata, key="sample")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=n_pcs)
else:
    sc.pp.neighbors(adata, n_pcs=n_pcs)
```

---

### 2.4 `rnaseq` — INCOMPLETE

**Status: No Snakefile exists.** Only `config_template.yaml`, `README.md`, and `samplesheet_template.csv` are present.

The config references these steps: FastQC → Trimmomatic → STAR (2-pass) → featureCounts → DESeq2 → clusterProfiler. None of these are implemented.

**What needs to be built:**
- `Snakefile` implementing all 6 steps
- `envs/` directory with conda environment YAML files for each tool group
- `scripts/` directory for DESeq2 and clusterProfiler R scripts

**Config audit (for when Snakefile is written):**

```yaml
# config_template.yaml — current settings:
threads_align: 16   # should be 64 for this hardware
mem_gb: 64          # STAR 2-pass for mammalian genome needs ~33 GB; 64 is fine
                    # but for large genomes (>5 GB) raise to 128
```

The config uses `threads_align` instead of the `threads` convention used by all other workflows. Standardize to `threads`.

**Missing from config:** No `threads_count:`, `threads_deseq2:`, or similar per-step thread controls. featureCounts is multi-threaded (`-T`); DESeq2 benefits from `BiocParallel`.

---

### 2.5 `atacseq` — INCOMPLETE

**Status: No Snakefile or config_template.yaml exists.** Only `samplesheet_template.csv` (identical format to rnaseq) is present.

The CLAUDE.md description lists: Bowtie2, MACS3, DiffBind. None are implemented.

**What needs to be built:**
- `config_template.yaml`
- `Snakefile` implementing: FastQC → Trimmomatic/fastp → Bowtie2 → samtools markdup → MACS3 → DiffBind → IDR
- `envs/` conda environments

**Hardware notes for when this is implemented:**
- Bowtie2 should use `threads: 32+`; `--no-discordant --no-mixed` for ATAC
- MACS3 is mostly single-threaded; parallelism via sample-level scatter
- DiffBind's DESeq2 backend benefits from BiocParallel

---

### 2.6 `genome_annotation` — INCOMPLETE

**Status: No Snakefile or config_template.yaml exists.** Only `samplesheet_template.csv` is present.

The CLAUDE.md description lists: BUSCO, RepeatModeler, BRAKER4, DIAMOND, InterPro. None are implemented.

**What needs to be built:**
- `config_template.yaml`
- `Snakefile` implementing the full annotation pipeline
- `envs/` conda environments

**Hardware notes for when this is implemented (highest priority workflow to implement):**

- **RepeatModeler** is the single longest-running step (~72–168 h for a 1 GB genome single-threaded). It supports `-pa` (number of parallel search jobs, each using multiple cores). For 160 cores: `-pa 20` with each job using 8 cores = 160 cores total. This is critical to configure correctly.
- **BRAKER4** runs AUGUSTUS + GeneMark in parallel. Should be given `--cores 80` minimum.
- **DIAMOND** nr database search: use `--threads 64 --block-size 10` for the 160-core machine (block-size 10 = ~100 GB RAM per thread block — check against 2.2 TB).
- **InterProScan** distributes jobs internally; pass `-cpu 64`.
- **BUSCO** uses Augustus internally which is single-threaded per gene, but BUSCO itself spawns parallel Augustus jobs: `--cpu 64`.
- The NR DIAMOND DB (`.../databases/nr.dmnd`) and InterProScan data already exist on disk per CLAUDE.md — paths should be referenced in config.

---

## 3. Priority Fix List (Ranked by Impact)

| Rank | Fix | Workflow | Impact | Effort |
|---|---|---|---|---|
| 1 | **Fix `PHENOTYPE` NameError** — replace `PHENOTYPE` with `METADATA` in `rule plots` input | `snp_association` | Crash at parse time — workflow is currently broken | 1 line |
| 2 | **Fix R input name mismatch** — `snakemake@input[["gwas"]]` → `snakemake@input[["gwas_salt"]]` in `manhattan_plot.R` | `snp_association` | R script crashes; plots never generated | 1 line |
| 3 | **Raise `threads` default to 64 / `mem_gb` to 256** in `wgs_snp` config | `wgs_snp` | BWA-MEM2 runs 4× faster; DeepVariant does not OOM at high shards | Config change |
| 4 | **Mark `markdup.bam` as `temp()`** in `wgs_snp` | `wgs_snp` | Saves 20–100 GB × N samples of stranded disk space | 2 lines |
| 5 | **Mark STARsolo BAM as `temp()`** in `scrnaseq` | `scrnaseq` | Saves 5–30 GB × N samples of disk space | 1 line |
| 6 | **Add `--threads` to all PLINK2 calls** and `resources: mem_mb` to all rules | `snp_association` | 4–8× speedup on VCF-to-BED, LD prune, PCA | ~10 lines |
| 7 | **Conditionally define `gwas_bodysize`** with `if HAS_COV:` | `snp_association` | Prevents running GEMMA on empty phenotype file | 2 lines |
| 8 | **Switch scDblFinder to `MulticoreParam`** + add `threads: 8` to rule | `scrnaseq` | 4–8× speedup on doublet detection for large samples | 3 lines |
| 9 | **Add `-Xmx16g` to snpEff ann** and add `--mem-gbytes` to glnexus_cli | `wgs_snp` | Prevents JVM OOM / GC thrash on large genomes; enforces GLnexus memory cap | 2 lines |
| 10 | **Add GPU resource sentinel** (`resources: gpu=1`) to DeepVariant rule and document `--resources gpu=1` at launch | `wgs_snp` | Prevents accidental parallel GPU jobs that would OOM VRAM | 2 lines + docs |

---

## 4. Estimated Runtime Table

Runtimes assume recommended settings applied (threads=64, mem_gb=256 for wgs_snp) and the hardware described above. All estimates assume fast NVMe I/O.

### `wgs_snp` (30 samples × 30× WGS, 1 GB genome)

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| fastp (per sample) | 16 threads × 30 parallel | ~20 min total |
| BWA-MEM2 align (per sample) | 64 threads, samples in batches of 2 | ~2–3 h total |
| samtools sort + markdup | 8 threads × parallel | ~30 min total |
| DeepVariant GPU (per sample, sequential) | 64 shards + GPU | ~25–40 min/sample → **12–20 h total** |
| GLnexus joint call | 64 threads | ~30–60 min |
| Filter + SNPeff | single-threaded mostly | ~15 min |
| **Total** | | **~14–24 h** |

DeepVariant is the dominant bottleneck. Running 2 samples concurrently with VRAM monitoring could halve this to ~7–12 h.

### `snp_association` (60 samples, 5M SNPs)

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| prepare_inputs | single-threaded | <5 min |
| VCF → PLINK (with --threads 16) | 16 threads | ~10–20 min |
| LD prune + PCA | 16 threads | ~5–10 min |
| Kinship matrix (GEMMA) | single-threaded | ~30–60 min |
| GWAS salt + bodysize (GEMMA) | single-threaded each | ~1–2 h each |
| Fst scan (VCFtools) | single-threaded | ~5–10 min |
| Plots + annotation | single-threaded | ~10–15 min |
| **Total** | | **~4–6 h** |

### `scrnaseq` (8 samples × 5,000 cells each)

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| FastQC (per sample) | 4 threads × 8 parallel | ~10 min |
| STAR index build | 32 threads | ~20–40 min |
| STARsolo align (per sample) | 32 threads × 2 parallel | ~20–30 min/sample → ~2–3 h total |
| scDblFinder (with MulticoreParam 8) | 8 threads × parallel | ~5–10 min/sample |
| Scanpy pipeline (all samples merged) | 32 threads (limited) | ~30–60 min |
| **Total** | | **~3–5 h** |

With GPU-accelerated Scanpy (rapids-singlecell) the final pipeline step drops to ~5–10 min.

### `rnaseq` (12 samples, paired-end, mammalian)

*(Estimate based on planned tools; no Snakefile exists yet)*

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| FastQC + Trimmomatic | 8–16 threads × parallel | ~20–30 min |
| STAR 2-pass align (per sample) | 64 threads × 2 parallel | ~20–30 min/sample → ~2–3 h total |
| featureCounts | 16 threads | ~10–20 min |
| DESeq2 | BiocParallel 16 | ~5–15 min |
| clusterProfiler | single-threaded | ~10–20 min |
| **Total** | | **~3–5 h** |

### `atacseq` (8 samples, 4 per group)

*(Estimate based on planned tools; no Snakefile exists yet)*

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| QC + trimming | parallel | ~15–20 min |
| Bowtie2 align (per sample) | 32 threads × parallel | ~20–40 min/sample → ~1–2 h total |
| samtools markdup | parallel | ~20 min |
| MACS3 peak call | per-sample parallel | ~5–10 min/sample |
| DiffBind / IDR | parallel | ~30–60 min |
| **Total** | | **~2–4 h** |

### `genome_annotation` (1 GB genome, vertebrate, with RNAseq hints)

*(Estimate based on planned tools; no Snakefile exists yet)*

| Step | Parallelism | Estimated Wall-Clock |
|---|---|---|
| BUSCO | 64 cores | ~2–4 h |
| RepeatModeler | `-pa 20` (160 cores) | **~48–96 h** (dominant step) |
| RepeatMasker | 32 threads | ~6–12 h |
| BRAKER4 | 80 cores | ~12–24 h |
| DIAMOND nr search | 64 threads | ~4–8 h |
| InterProScan | 64 CPUs | ~8–16 h |
| **Total (serial)** | | **~80–160 h** |
| **Total (parallel where possible)** | | **~50–100 h** |

RepeatModeler is the critical path. Correct `-pa` configuration is essential.

---

## 5. Summary of All Code Changes Required

### Immediate (bugs that cause crashes):

**File:** `/home/cylin/hermes/workflows/snp_association/Snakefile` line 281
```python
# BEFORE:
pheno = PHENOTYPE,
# AFTER:
pheno = METADATA,
```

**File:** `/home/cylin/hermes/workflows/snp_association/scripts/manhattan_plot.R` line 9
```r
# BEFORE:
gwas_file   <- snakemake@input[["gwas"]]
# AFTER:
gwas_file   <- snakemake@input[["gwas_salt"]]
```

### High-priority performance fixes:

**File:** `/home/cylin/hermes/workflows/wgs_snp/config_template.yaml`
```yaml
# BEFORE:
threads: 16
mem_gb: 64
# AFTER:
threads: 64
mem_gb: 256
```

**File:** `/home/cylin/hermes/workflows/wgs_snp/Snakefile` lines 139-140
```python
# BEFORE:
bam = str(OUTDIR / "alignment" / "{sample}.markdup.bam"),
bai = str(OUTDIR / "alignment" / "{sample}.markdup.bam.bai"),
# AFTER:
bam = temp(str(OUTDIR / "alignment" / "{sample}.markdup.bam")),
bai = temp(str(OUTDIR / "alignment" / "{sample}.markdup.bam.bai")),
```

**File:** `/home/cylin/hermes/workflows/scrnaseq/Snakefile` line 127
```python
# BEFORE:
bam = str(OUTDIR / "starsolo" / "{sample}" / "Aligned.sortedByCoord.out.bam"),
# AFTER:
bam = temp(str(OUTDIR / "starsolo" / "{sample}" / "Aligned.sortedByCoord.out.bam")),
```

**File:** `/home/cylin/hermes/workflows/scrnaseq/scripts/doublet_detection.R` line 18
```r
# BEFORE:
sce <- scDblFinder(sce, BPPARAM = SerialParam())
# AFTER:
nworkers <- if (!is.null(snakemake@threads)) as.integer(snakemake@threads) else 4L
sce <- scDblFinder(sce, BPPARAM = MulticoreParam(workers = nworkers))
```

**File:** `/home/cylin/hermes/workflows/snp_association/Snakefile` — wrap `gwas_bodysize` in conditional:
```python
# BEFORE:
rule gwas_bodysize:
    ...

# AFTER:
if HAS_COV:
    rule gwas_bodysize:
        ...
```

---

*End of audit. Three workflows (rnaseq, atacseq, genome_annotation) require Snakefile implementation before they can be used. Two workflows (wgs_snp, scrnaseq) are functional but contain disk-waste and sub-optimal thread allocation. One workflow (snp_association) has two runtime-crash bugs that must be fixed before first use.*
