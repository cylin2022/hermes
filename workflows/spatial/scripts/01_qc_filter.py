"""
01_qc_filter.py — Load spatial data, compute QC metrics, filter spots/cells.

Supports: 10x Visium, Visium HD, Xenium
Output  : QC-filtered AnnData (.h5ad) + QC plots
"""

import sys
import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

data_dir     = sys.argv[1]
platform     = sys.argv[2]    # visium | visium_hd | xenium
sample_name  = sys.argv[3]
outdir       = Path(sys.argv[4])
min_counts   = int(sys.argv[5])
min_genes    = int(sys.argv[6])
max_counts   = int(sys.argv[7])
max_mt_pct   = float(sys.argv[8])
mt_prefix    = sys.argv[9]

sc.settings.verbosity = 2
outdir.mkdir(parents=True, exist_ok=True)
plot_dir = outdir / "plots"
plot_dir.mkdir(exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print(f"[{sample_name}] Loading {platform} data from {data_dir}...")

if platform in ("visium", "visium_hd"):
    adata = sc.read_visium(data_dir)
    adata.var_names_make_unique()
    adata.obs_names = [f"{sample_name}_{b}" for b in adata.obs_names]
elif platform == "xenium":
    # Xenium In Situ output: cell_feature_matrix/ directory
    try:
        import spatialdata_io
    except ImportError:
        raise ImportError(
            "spatialdata-io is required for Xenium data. "
            "Install with: pip install spatialdata-io"
        )
    sdata = spatialdata_io.xenium(data_dir)
    # spatialdata stores tables under sdata.tables; key may vary by version
    table_key = "table" if "table" in sdata.tables else list(sdata.tables.keys())[0]
    adata = sdata.tables[table_key].copy()
    adata.obs_names = [f"{sample_name}_{c}" for c in adata.obs_names]
else:
    raise ValueError(f"Unknown platform: {platform}. Use visium | visium_hd | xenium")

adata.obs["sample"] = sample_name
print(f"  Raw: {adata.n_obs} spots × {adata.n_vars} genes")

# ── QC metrics ────────────────────────────────────────────────────────────────
adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# ── Violin / spatial QC plots ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, metric, label in zip(
    axes,
    ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    ["Total UMI", "Genes detected", "MT%"],
):
    ax.violinplot(adata.obs[metric], showmedians=True)
    ax.set_title(label)
    ax.set_xticks([])
plt.suptitle(f"{sample_name} — before filtering")
plt.tight_layout()
plt.savefig(plot_dir / "qc_violin_before.pdf", bbox_inches="tight")
plt.close()

# Spatial QC maps
if hasattr(adata, "obsm") and "spatial" in adata.obsm:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric in zip(axes,
                          ["total_counts", "n_genes_by_counts", "pct_counts_mt"]):
        sq.pl.spatial_scatter(adata, color=metric, ax=ax, show=False,
                               colorbar=True, title=metric)
    plt.suptitle(f"{sample_name} — QC spatial distribution")
    plt.tight_layout()
    plt.savefig(plot_dir / "qc_spatial.pdf", bbox_inches="tight")
    plt.close()

# ── Filter ────────────────────────────────────────────────────────────────────
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_cells(adata, max_counts=max_counts)
sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct].copy()
sc.pp.filter_genes(adata, min_cells=10)

n_after = adata.n_obs
print(f"  After QC: {n_after} spots retained ({n_before - n_after} removed, "
      f"{100 * n_after / n_before:.1f}% pass rate)")
print(f"  Genes after filter: {adata.n_vars}")

# ── Save raw counts (needed for downstream deconvolution) ─────────────────────
adata.layers["counts"] = adata.X.copy()

out_h5ad = outdir / f"{sample_name}_qc.h5ad"
adata.write_h5ad(out_h5ad)
print(f"  Saved: {out_h5ad}")

# QC summary text
with open(outdir / f"{sample_name}_qc_summary.txt", "w") as f:
    f.write(f"Sample:        {sample_name}\n")
    f.write(f"Platform:      {platform}\n")
    f.write(f"Spots before:  {n_before}\n")
    f.write(f"Spots after:   {n_after}\n")
    f.write(f"Genes after:   {adata.n_vars}\n")
    f.write(f"Median UMI:    {adata.obs['total_counts'].median():.0f}\n")
    f.write(f"Median genes:  {adata.obs['n_genes_by_counts'].median():.0f}\n")
    f.write(f"Median MT%:    {adata.obs['pct_counts_mt'].median():.2f}\n")
