"""
Scanpy pipeline: load STARsolo counts → QC → normalise → HVG → PCA → UMAP →
Leiden clustering → marker gene detection.
"""

import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# ── Parameters from Snakemake ────────────────────────────────────────────────
samples      = snakemake.params.samples
starsolo_dir = Path(snakemake.params.starsolo_dir)
doublet_dir  = Path(snakemake.params.doublet_dir)
min_genes    = snakemake.params.min_genes
max_genes    = snakemake.params.max_genes
max_mt_pct   = snakemake.params.max_mt_pct
mt_prefix    = snakemake.params.mt_prefix
n_hvg        = snakemake.params.n_hvg
n_pcs        = snakemake.params.n_pcs
resolution   = snakemake.params.resolution
outdir       = Path(snakemake.params.outdir)
n_jobs       = snakemake.threads

sc.settings.n_jobs = n_jobs
sc.settings.verbosity = 2
outdir.mkdir(parents=True, exist_ok=True)

# ── Load each sample ─────────────────────────────────────────────────────────
adatas = []
for sample in samples:
    mtx_dir = starsolo_dir / sample / "Solo.out" / "Gene" / "filtered"
    print(f"[scanpy] Loading {sample} from {mtx_dir}")
    adata = sc.read_10x_mtx(str(mtx_dir), var_names="gene_ids", cache=True)
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]
    adata.obs["sample"] = sample

    doublet_file = doublet_dir / f"{sample}_doublet_scores.csv"
    if doublet_file.exists():
        db = pd.read_csv(doublet_file, index_col="barcode")
        adata.obs = adata.obs.join(db, how="left")
    adatas.append(adata)

if len(adatas) == 1:
    adata = adatas[0]
else:
    adata = sc.concat(adatas, label="sample", keys=samples, join="outer", fill_value=0)

print(f"[scanpy] Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# ── Remove doublets ───────────────────────────────────────────────────────────
if "doublet_class" in adata.obs.columns:
    n_before = adata.n_obs
    adata = adata[adata.obs["doublet_class"] != "doublet"].copy()
    print(f"[scanpy] Removed {n_before - adata.n_obs} doublets → {adata.n_obs} cells remain")

# ── QC metrics ───────────────────────────────────────────────────────────────
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# QC violin plot
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, key in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
    ax.violinplot(adata.obs[key].values, showmedians=True)
    ax.set_title(key)
fig.savefig(str(outdir / "qc_violin.pdf"), bbox_inches="tight")
plt.close()

# ── Cell filtering ────────────────────────────────────────────────────────────
sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs.n_genes_by_counts < max_genes].copy()
adata = adata[adata.obs.pct_counts_mt < max_mt_pct].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"[scanpy] After QC filter: {adata.n_obs} cells × {adata.n_vars} genes")

# ── Normalisation + log ───────────────────────────────────────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep raw counts for DE

# ── HVG → PCA → neighbours → UMAP → Leiden ───────────────────────────────────
batch_key = "sample" if len(samples) > 1 else None
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=batch_key)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
sc.pp.neighbors(adata, n_pcs=n_pcs)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=resolution)
print(f"[scanpy] Leiden clusters: {adata.obs['leiden'].nunique()} (resolution={resolution})")

# ── Marker genes ──────────────────────────────────────────────────────────────
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True, n_genes=100)
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(str(outdir / "marker_genes.csv"), index=False)

# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata, color="sample",  ax=axes[1], show=False, title="Sample")
fig.savefig(str(outdir / "umap_clusters.pdf"), bbox_inches="tight")
plt.close()

# Top 3 markers per cluster dot plot
top_markers = (
    markers.groupby("group")
    .apply(lambda x: x.nlargest(3, "scores")["names"].tolist())
    .explode()
    .unique()
    .tolist()
)
# Filter to genes present in adata
top_markers = [g for g in top_markers if g in adata.raw.var_names][:50]
if top_markers:
    fig = sc.pl.dotplot(adata, top_markers, groupby="leiden", use_raw=True,
                        show=False, return_fig=True)
    fig.savefig(str(outdir / "dotplot_markers.pdf"), bbox_inches="tight")
    plt.close()

# ── Save ──────────────────────────────────────────────────────────────────────
adata.write(str(outdir / "final_adata.h5ad"), compression="gzip")
print(f"[scanpy] Done. Saved to {outdir / 'final_adata.h5ad'}")
