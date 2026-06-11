# scRNA-seq Analysis Workflow

Single-cell RNA-seq: raw reads → QC → cell barcode calling → doublet removal → clustering → marker genes

## Pipeline steps

1. **FastQC + MultiQC** — raw read quality assessment
2. **STARsolo** — genome alignment + cell barcode/UMI demultiplexing (open-source CellRanger equivalent)
3. **scDblFinder** (R/Bioconductor) — doublet detection per sample
4. **Scanpy pipeline**
   - Merge all samples into one AnnData object
   - QC filtering (min/max genes per cell, % mitochondrial reads)
   - Remove doublets
   - Normalize (CPM) + log1p
   - Highly variable gene selection
   - PCA → KNN graph → UMAP
   - Leiden clustering
   - Wilcoxon rank-sum marker genes per cluster

## Supported chemistries

| Value | Platform |
|-------|----------|
| `10xv3` | 10x Chromium v3 (default) |
| `10xv2` | 10x Chromium v2 |
| `dropseq` | Drop-seq |

## Key outputs

| File | Description |
|------|-------------|
| `scanpy/final_adata.h5ad` | Fully processed AnnData (load in Scanpy/Seurat) |
| `scanpy/umap_clusters.pdf` | UMAP coloured by cluster and sample |
| `scanpy/marker_genes.csv` | Top 100 marker genes per cluster (Wilcoxon) |
| `scanpy/dotplot_markers.pdf` | Dot plot of top 3 markers per cluster |
| `scanpy/qc_violin.pdf` | QC violin plots (genes, counts, % MT) |
| `doublets/{sample}_doublet_scores.csv` | Per-cell doublet scores |
| `multiqc/multiqc_report.html` | Combined QC report |

## Notes for non-model organisms

- Set `mt_gene_prefix` to match your genome annotation (common: `MT-`, `mt:`, `Mt`, `chrM_`)
- Set `star_sa_index_nbases` based on genome size:
  - Human/mouse (~3 Gb): 14
  - Fish/frog (~1 Gb): 13
  - Insect (~300 Mb): 12–13
  - Small genomes (<50 Mb): 10–11

## Estimated runtime (160 cores, 4 samples, human)

| Step | Time |
|------|------|
| STAR index build | ~1 h |
| STARsolo (per sample, 50k cells) | ~30 min |
| Doublet detection | ~10 min/sample |
| Scanpy pipeline | ~30 min |
| **Total** | **~3–4 h** |
