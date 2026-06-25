"""
04_deconvolution.py — Cell type deconvolution with cell2location.

cell2location uses a Bayesian model to infer the number of cells of each type
per Visium spot, given a reference scRNA-seq signature.

Input:
  - Spatial AnnData with raw counts in .layers["counts"]
  - Reference scRNA-seq AnnData with cell type labels in .obs[celltype_col]
Output:
  - AnnData with cell type abundance in .obs (q05_cell_abundance_w_sf_*)
  - Spatial maps of cell type abundance
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

h5ad_spatial   = sys.argv[1]
h5ad_reference = sys.argv[2]
celltype_col   = sys.argv[3]
sample_name    = sys.argv[4]
outdir         = Path(sys.argv[5])
n_cells_per_loc = int(sys.argv[6])     # expected cells per spot
detection_alpha = float(sys.argv[7])
n_threads       = int(sys.argv[8])

import torch
import cell2location
from cell2location.models import RegressionModel, Cell2location

sc.settings.verbosity = 2
plot_dir = outdir / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

# cell2location uses scvi-tools internally; newer versions use accelerator= not use_gpu=
gpu_available = torch.cuda.is_available()
print(f"[{sample_name}] GPU available: {gpu_available}")

import scvi
import re as _re
# Handle pre-release version strings like "0.20.0a1" or "1.0.0.post1"
_scvi_ver_match = _re.match(r"(\d+)\.(\d+)", scvi.__version__)
if _scvi_ver_match:
    _scvi_version = (int(_scvi_ver_match.group(1)), int(_scvi_ver_match.group(2)))
else:
    _scvi_version = (0, 0)
# scvi-tools ≥ 0.20 uses accelerator/devices instead of use_gpu=
if _scvi_version >= (0, 20):
    _train_kwargs = {"accelerator": "gpu" if gpu_available else "cpu", "devices": 1}
else:
    _train_kwargs = {"use_gpu": gpu_available}

# ── Load data ─────────────────────────────────────────────────────────────────
adata_vis = sc.read_h5ad(h5ad_spatial)
adata_ref = sc.read_h5ad(h5ad_reference)

# Preflight: reference must be non-empty and have the required cell type column
if adata_ref.n_obs == 0:
    raise ValueError(
        f"Reference h5ad '{h5ad_reference}' has 0 cells. "
        "Provide a non-empty scRNA-seq reference."
    )
if celltype_col not in adata_ref.obs.columns:
    raise ValueError(
        f"Cell type column '{celltype_col}' not found in reference .obs. "
        f"Available columns: {list(adata_ref.obs.columns)}"
    )

# Preflight: reference counts layer check — cell2location requires raw integer counts
if "counts" not in adata_ref.layers:
    import scipy.sparse as _sp
    X = adata_ref.X
    X_dense = X.toarray() if _sp.issparse(X) else X
    if not np.allclose(X_dense, np.round(X_dense)):
        raise ValueError(
            "Reference adata.X does not appear to contain raw integer counts "
            "and no 'counts' layer was found. "
            "Please provide raw (un-normalized) counts in .X or .layers['counts']."
        )
    print("  Warning: no 'counts' layer in reference; using .X (appears to be integer counts)")
else:
    print("  Found 'counts' layer in reference; will use it for cell2location.")

print(f"  Spatial: {adata_vis.n_obs} spots × {adata_vis.n_vars} genes")
print(f"  Reference: {adata_ref.n_obs} cells × {adata_ref.n_vars} genes")
print(f"  Cell types: {adata_ref.obs[celltype_col].nunique()}")

# ── Step 1: Estimate reference signatures ─────────────────────────────────────
# Filter reference to genes present in spatial data
shared_genes = adata_vis.var_names.intersection(adata_ref.var_names)
adata_ref_sub = adata_ref[:, shared_genes].copy()
adata_vis_sub = adata_vis[:, shared_genes].copy()

print(f"  Shared genes: {len(shared_genes)}")

# cell2location requires raw integer counts
if "counts" in adata_ref.layers:
    adata_ref_sub.X = adata_ref_sub.layers["counts"]

# Filter lowly expressed genes in reference
sc.pp.filter_genes(adata_ref_sub, min_cells=10)

# Prepare for regression model
cell2location.models.RegressionModel.setup_anndata(
    adata_ref_sub,
    labels_key = celltype_col,
    batch_key  = None,
)
mod_ref = RegressionModel(adata_ref_sub)
mod_ref.train(max_epochs=250, batch_size=2500, **_train_kwargs)

# Extract signatures
adata_ref_sub = mod_ref.export_posterior(
    adata_ref_sub,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500, **_train_kwargs},
)
inf_aver = adata_ref_sub.varm["means_per_cluster_mu_fg"][
    [f"means_per_cluster_mu_fg_{c}" for c in adata_ref_sub.uns["mod"]["factor_names"]]
].copy()
inf_aver.columns = adata_ref_sub.uns["mod"]["factor_names"]

print(f"  Reference signature estimated for {inf_aver.shape[1]} cell types")

# ── Step 2: Map to spatial data ───────────────────────────────────────────────
if "counts" in adata_vis_sub.layers:
    adata_vis_sub.X = adata_vis_sub.layers["counts"]

# Align genes
intersect = np.intersect1d(adata_vis_sub.var_names, inf_aver.index)
adata_vis_sub = adata_vis_sub[:, intersect].copy()
inf_aver_sub  = inf_aver.loc[intersect]

cell2location.models.Cell2location.setup_anndata(
    adata_vis_sub, batch_key=None
)
mod_c2l = Cell2location(
    adata_vis_sub,
    cell_state_df           = inf_aver_sub,
    N_cells_per_location    = n_cells_per_loc,
    detection_alpha         = detection_alpha,
)
mod_c2l.train(max_epochs=30000, batch_size=2500, train_size=1, **_train_kwargs)

adata_vis_sub = mod_c2l.export_posterior(
    adata_vis_sub,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500, **_train_kwargs},
)

# ── Extract and store cell type abundances ────────────────────────────────────
# Use 5th percentile estimates (conservative, avoids noise)
ctypes = inf_aver_sub.columns.tolist()
# Assert expected output key exists after export_posterior
assert "q05_cell_abundance_w_sf" in adata_vis_sub.obsm, (
    "cell2location export_posterior did not produce 'q05_cell_abundance_w_sf' in .obsm. "
    f"Available keys: {list(adata_vis_sub.obsm.keys())}"
)
if "q05_cell_abundance_w_sf" in adata_vis_sub.obsm:
    cell_abund = adata_vis_sub.obsm["q05_cell_abundance_w_sf"]
    # cell2location column names may be prefixed; force rename to clean cell type names
    if len(cell_abund.columns) == len(ctypes):
        cell_abund = cell_abund.copy()
        cell_abund.columns = ctypes
    adata_vis.obs[ctypes] = cell_abund.values

# ── Save model ────────────────────────────────────────────────────────────────
mod_c2l.save(str(outdir / "cell2location_model"), overwrite=True)

out_h5ad = outdir / f"{sample_name}_deconv.h5ad"
adata_vis.write_h5ad(out_h5ad)
print(f"  Saved: {out_h5ad}")

# ── Spatial abundance maps (top 8 cell types by mean abundance) ───────────────
if "spatial" in adata_vis.obsm and ctypes:
    top_ct = (adata_vis.obs[ctypes].mean().sort_values(ascending=False)
              .head(8).index.tolist())
    fig = sq.pl.spatial_scatter(
        adata_vis,
        color      = top_ct,
        ncols      = 4,
        cmap       = "Reds",
        show       = False,
        title      = top_ct,
        return_fig = True,
    )
    fig.suptitle(f"{sample_name} — Cell type abundance (cell2location)", y=1.02)
    fig.savefig(plot_dir / "spatial_celltype_abundance.pdf", bbox_inches="tight")
    plt.close(fig)

    # Stacked bar per cluster
    if "leiden" in adata_vis.obs.columns:
        cluster_ct = adata_vis.obs.groupby("leiden")[ctypes].mean()
        cluster_ct_norm = cluster_ct.div(cluster_ct.sum(axis=1), axis=0)
        cluster_ct_norm.plot(kind="bar", stacked=True, figsize=(12, 5),
                              colormap="tab20")
        plt.ylabel("Proportion")
        plt.title(f"{sample_name} — Cell type composition per cluster")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        plt.savefig(plot_dir / "cluster_celltype_composition.pdf",
                    bbox_inches="tight")
        plt.close()

print(f"  Deconvolution complete. {len(ctypes)} cell types mapped.")
