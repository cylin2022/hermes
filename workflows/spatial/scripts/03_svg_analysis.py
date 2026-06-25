"""
03_svg_analysis.py — Spatially Variable Genes via Moran's I (Squidpy).

Moran's I: global spatial autocorrelation statistic.
  I ≈ +1 → spatially clustered expression (SVG)
  I ≈  0 → random spatial distribution
  I ≈ -1 → spatially dispersed (checkerboard)

Output: SVG table (CSV) + spatial expression maps of top genes.
"""

import sys
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

h5ad_in    = sys.argv[1]
sample_name = sys.argv[2]
outdir     = Path(sys.argv[3])
n_svg      = int(sys.argv[4])     # top N SVGs to report
n_jobs     = int(sys.argv[5])     # parallel jobs

sc.settings.verbosity = 2
plot_dir = outdir / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(h5ad_in)
print(f"[{sample_name}] Computing Moran's I for {adata.n_vars} genes...")

# ── Rebuild spatial graph if not present ──────────────────────────────────────
if "spatial_connectivities" not in adata.obsp:
    sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")

# ── Moran's I ─────────────────────────────────────────────────────────────────
# Use log-normalized layer; limit to HVGs for speed if dataset is large
if "highly_variable" in adata.var.columns:
    genes_to_test = adata.var_names[adata.var["highly_variable"]].tolist()
else:
    genes_to_test = adata.var_names.tolist()

sq.gr.spatial_autocorr(
    adata,
    mode      = "moran",
    genes     = genes_to_test,
    n_perms   = 1000,       # permutation test for p-value
    corr_method = "fdr_bh", # Benjamini-Hochberg FDR correction
    n_jobs    = n_jobs,
    layer     = "log1p",
)

# ── Extract results ───────────────────────────────────────────────────────────
svg_df = adata.uns["moranI"].copy()
svg_df = svg_df.sort_values("I", ascending=False)
svg_df.index.name = "gene"

print(f"  Total genes tested: {len(svg_df)}")
print(f"  FDR < 0.05: {(svg_df['pval_norm_fdr_bh'] < 0.05).sum()}")
print(f"  Top 5 SVGs:")
print(svg_df.head(5)[["I", "pval_norm", "pval_norm_fdr_bh"]].to_string())

top_svgs = svg_df[svg_df["pval_norm_fdr_bh"] < 0.05].head(n_svg)
if len(top_svgs) == 0:
    # Fallback: no gene passed FDR threshold; report top by raw Moran's I score
    print("  WARNING: No genes passed FDR < 0.05. Reporting top genes by I score.")
    top_svgs = svg_df.head(n_svg)
svg_df.to_csv(outdir / f"{sample_name}_moranI_all.csv")
top_svgs.to_csv(outdir / f"{sample_name}_svg_top.csv")

# ── Spatial expression maps of top 12 SVGs ────────────────────────────────────
top12 = top_svgs.head(12).index.tolist()
if top12 and "spatial" in adata.obsm:
    fig = sq.pl.spatial_scatter(
        adata,
        color      = top12,
        layer      = "log1p",
        ncols      = 4,
        show       = False,
        title      = [f"{g}\nI={svg_df.loc[g,'I']:.3f}" for g in top12],
        return_fig = True,
    )
    fig.suptitle(f"{sample_name} — Top SVGs (Moran's I)", fontsize=12, y=1.02)
    fig.savefig(plot_dir / "spatial_top_svgs.pdf", bbox_inches="tight")
    plt.close(fig)

# ── Moran's I score distribution ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(svg_df["I"], bins=50, color="#4E79A7", edgecolor="white", linewidth=0.5)
ax.axvline(0, color="grey", linestyle="--", linewidth=1)
ax.set_xlabel("Moran's I")
ax.set_ylabel("Number of genes")
ax.set_title(f"{sample_name} — Moran's I distribution")
plt.tight_layout()
plt.savefig(plot_dir / "moranI_distribution.pdf", bbox_inches="tight")
plt.close()

print(f"  SVG analysis complete. Results: {outdir / f'{sample_name}_svg_top.csv'}")
