"""
05_spatial_stats.py — Spatial statistics with Squidpy.

Analyses:
  1. Neighborhood enrichment  — which cluster pairs co-localize on tissue?
  2. Co-occurrence score      — distance-dependent cluster co-occurrence
  3. Ligand-receptor (LIANA)  — cell communication from spatial context
  4. Ripley's statistics      — point pattern analysis (optional)

Output: stat tables + heatmaps + interaction plots
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

h5ad_in               = sys.argv[1]
sample_name           = sys.argv[2]
outdir                = Path(sys.argv[3])
run_nhood_enrichment  = sys.argv[4].lower() == "true"
run_co_occurrence     = sys.argv[5].lower() == "true"
run_ripley            = sys.argv[6].lower() == "true"
co_occ_interval       = int(sys.argv[7])
n_threads             = int(sys.argv[8])

sc.settings.verbosity = 2
plot_dir = outdir / "plots"
plot_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(h5ad_in)
print(f"[{sample_name}] Spatial statistics — {adata.n_obs} spots, "
      f"{adata.obs['leiden'].nunique()} clusters")

# Ensure spatial graph exists
if "spatial_connectivities" not in adata.obsp:
    sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")

# ── 1. Neighborhood enrichment ────────────────────────────────────────────────
if run_nhood_enrichment:
    print("  Running neighborhood enrichment...")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden", seed=42, n_perms=1000)

    fig, ax = plt.subplots(figsize=(8, 6))
    sq.pl.nhood_enrichment(
        adata, cluster_key="leiden",
        method="zscore",
        title=f"{sample_name} — Neighborhood enrichment",
        ax=ax, show=False,
    )
    plt.tight_layout()
    plt.savefig(plot_dir / "nhood_enrichment.pdf", bbox_inches="tight")
    plt.close()

    # Save z-score matrix
    zscore = adata.uns["leiden_nhood_enrichment"]["zscore"]
    clusters = adata.obs["leiden"].cat.categories.tolist()
    pd.DataFrame(zscore, index=clusters, columns=clusters).to_csv(
        outdir / f"{sample_name}_nhood_enrichment.csv"
    )
    print("    Done.")

# ── 2. Co-occurrence ─────────────────────────────────────────────────────────
if run_co_occurrence:
    print("  Running co-occurrence analysis...")
    sq.gr.co_occurrence(
        adata,
        cluster_key = "leiden",
        interval    = co_occ_interval,
        n_jobs      = n_threads,
        seed        = 42,
    )
    fig = sq.pl.co_occurrence(
        adata,
        cluster_key = "leiden",
        clusters    = adata.obs["leiden"].cat.categories[:6].tolist(),
        show        = False,
        return_fig  = True,
    )
    fig.suptitle(f"{sample_name} — Co-occurrence score", y=1.01)
    fig.savefig(plot_dir / "co_occurrence.pdf", bbox_inches="tight")
    plt.close(fig)
    print("    Done.")

# ── 3. Centrality scores ──────────────────────────────────────────────────────
print("  Computing centrality scores...")
sq.gr.centrality_scores(adata, cluster_key="leiden")
fig = sq.pl.centrality_scores(
    adata, cluster_key="leiden",
    show=False,
    return_fig=True,
)
fig.suptitle(f"{sample_name} — Centrality scores", y=1.01)
fig.savefig(plot_dir / "centrality_scores.pdf", bbox_inches="tight")
plt.close(fig)

# ── 4. Ripley statistics (optional, slow) ────────────────────────────────────
if run_ripley:
    print("  Computing Ripley's K/L statistics...")
    sq.gr.ripley(adata, cluster_key="leiden", mode="L", n_simulations=100,
                 n_jobs=n_threads)
    fig = sq.pl.ripley(adata, cluster_key="leiden", mode="L", show=False,
                       return_fig=True)
    fig.suptitle(f"{sample_name} — Ripley's L")
    fig.savefig(plot_dir / "ripley_L.pdf", bbox_inches="tight")
    plt.close(fig)

# ── 5. Ligand-receptor with LIANA ─────────────────────────────────────────────
# LIANA ≥1.0 uses liana.mt.rank_aggregate(); older versions used by_sample().
# We try the modern API first and fall back gracefully.
try:
    import liana
    print("  Running LIANA ligand-receptor analysis...")

    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"]

    liana.mt.rank_aggregate(
        adata,
        groupby       = "leiden",
        resource_name = "consensus",
        use_raw       = False,
        verbose       = True,
        n_perms       = 100,
    )

    # liana result is stored in adata.uns["liana_res"] as a DataFrame
    if "liana_res" in adata.uns:
        liana_df = adata.uns["liana_res"]
        if hasattr(liana_df, "sort_values"):
            liana_df.sort_values("aggregate_rank").head(500).to_csv(
                outdir / f"{sample_name}_liana_top500.csv", index=False
            )

        # Dot plot of top interactions using LIANA's built-in plotter
        try:
            fig = liana.pl.dotplot(
                adata,
                colour      = "lr_means",
                size        = "cellphone_pvals",
                top_n       = 30,
                return_fig  = True,
            )
            fig.savefig(plot_dir / "liana_dotplot.pdf", bbox_inches="tight")
            plt.close("all")
        except Exception as pe:
            print(f"    LIANA plot skipped ({pe})")

    print("    LIANA done.")

except ImportError:
    print("  LIANA not installed — skipping L-R analysis")
except Exception as e:
    print(f"  LIANA failed ({e}) — skipping")

# ── Save final AnnData ────────────────────────────────────────────────────────
out_h5ad = outdir / f"{sample_name}_final.h5ad"
adata.write_h5ad(out_h5ad)
print(f"  Saved final AnnData: {out_h5ad}")
print(f"  Spatial stats complete.")
