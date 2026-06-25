"""
02_normalize_cluster.py — Normalize, HVG, PCA, UMAP, Leiden clustering.
Output: clustered AnnData + UMAP + spatial cluster plots.
"""

import sys
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

h5ad_in     = sys.argv[1]
sample_name = sys.argv[2]
outdir      = Path(sys.argv[3])
n_hvg       = int(sys.argv[4])
n_pcs       = int(sys.argv[5])
n_neighbors = int(sys.argv[6])
resolution  = float(sys.argv[7])
spatial_k   = int(sys.argv[8])   # spatial neighbor count for Squidpy graph
seed        = int(sys.argv[9])   # random seed for reproducibility

sc.settings.verbosity = 2
plot_dir = outdir / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(h5ad_in)
print(f"[{sample_name}] Loaded {adata.n_obs} spots × {adata.n_vars} genes")

# ── Normalize ─────────────────────────────────────────────────────────────────
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

# ── HVG + PCA ─────────────────────────────────────────────────────────────────
sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3",
                             layer="counts")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True)

# ── kNN graph + UMAP ──────────────────────────────────────────────────────────
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
sc.tl.umap(adata, random_state=seed)

# ── Leiden clustering ─────────────────────────────────────────────────────────
sc.tl.leiden(adata, resolution=resolution, key_added="leiden", random_state=seed)
n_clusters = adata.obs["leiden"].nunique()
print(f"  Leiden clusters: {n_clusters} (resolution={resolution})")

# ── Build spatial neighbor graph (Squidpy) ─────────────────────────────────────
# coord_type="generic" uses kNN on pixel coordinates (correct for Visium + Xenium).
# coord_type="visium" is Visium-specific hex grid — but pixel coords work for both.
sq.gr.spatial_neighbors(adata, n_neighs=spatial_k, coord_type="generic")

# ── Marker genes per cluster ──────────────────────────────────────────────────
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                         layer="log1p", use_raw=False)

# ── Plots ─────────────────────────────────────────────────────────────────────
# Capture Scanpy figures explicitly to avoid plt.savefig() grabbing wrong figure.
fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(adata, color="leiden", ax=ax, show=False)
fig.savefig(plot_dir / "umap_leiden.pdf", bbox_inches="tight")
plt.close(fig)

# Spatial distribution of clusters
if "spatial" in adata.obsm:
    fig = sq.pl.spatial_scatter(adata, color="leiden", show=False,
                                 title=f"{sample_name} — Leiden clusters",
                                 return_fig=True)
    fig.savefig(plot_dir / "spatial_leiden.pdf", bbox_inches="tight")
    plt.close(fig)

    # Top QC metrics on tissue
    for metric in ["total_counts", "n_genes_by_counts"]:
        fig = sq.pl.spatial_scatter(adata, color=metric, show=False,
                                     cmap="viridis", title=metric,
                                     return_fig=True)
        fig.savefig(plot_dir / f"spatial_{metric}.pdf", bbox_inches="tight")
        plt.close(fig)

# Dot plot: top 5 markers per cluster
axes = sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False,
                                        groupby="leiden", standard_scale="var",
                                        return_fig=True)
axes.savefig(plot_dir / "dotplot_markers.pdf", bbox_inches="tight")
plt.close("all")

# ── Save ──────────────────────────────────────────────────────────────────────
out_h5ad = outdir / f"{sample_name}_clustered.h5ad"
adata.write_h5ad(out_h5ad)
print(f"  Saved: {out_h5ad}")

# Cluster stats
with open(outdir / f"{sample_name}_cluster_summary.txt", "w") as f:
    f.write(f"Clusters: {n_clusters}\n")
    f.write(adata.obs["leiden"].value_counts().to_string())
    f.write("\n")
