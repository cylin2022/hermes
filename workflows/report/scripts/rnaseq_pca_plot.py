"""PCA plot for rnaseq report.

The rnaseq workflow does not produce a vst_counts.tsv matrix directly.
Instead it writes {contrast}_rlog.rds files per contrast. Since reading R RDS
objects in Python requires rpy2 (not always available), this script produces a
placeholder figure noting that PCA requires re-extraction from the RDS objects.

To enable a real PCA plot, add an R script that exports the rlog matrix to TSV
and update this script accordingly.
"""
import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

input_files = snakemake.input   # list of *_rlog.rds paths (may be empty)
out         = snakemake.output[0]
deseq2_dir  = Path(snakemake.params.deseq2_dir)

# Attempt to find a TSV-format count matrix exported alongside the rlog RDS
vst_candidates = list(deseq2_dir.glob("*_vst_counts.tsv")) + list(deseq2_dir.glob("vst_counts.tsv"))

if vst_candidates:
    import pandas as pd
    import numpy as np
    from sklearn.decomposition import PCA

    df = pd.read_csv(vst_candidates[0], sep="\t", index_col=0)
    pca = PCA(n_components=min(2, df.shape[1]))
    coords = pca.fit_transform(df.T)
    var_exp = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(coords[:, 0], coords[:, 1] if coords.shape[1] > 1 else [0] * len(coords),
               s=80, alpha=0.8)
    for i, name in enumerate(df.columns):
        ax.annotate(name, (coords[i, 0], coords[i, 1] if coords.shape[1] > 1 else 0),
                    fontsize=8, ha="center", va="bottom")
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)" if len(var_exp) > 1 else "PC2")
    ax.set_title("PCA — VST normalised counts")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"PCA plot saved from {vst_candidates[0]}")
else:
    # No VST matrix available — produce informative placeholder
    warnings.warn(
        "PCA plot skipped: rnaseq workflow does not export a vst_counts.tsv matrix. "
        "Add an R export step or use MultiQC PCA instead."
    )
    rds_files = list(input_files) if input_files else []
    n_rds = len(rds_files)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.6, "PCA plot not available",
            ha="center", va="center", transform=ax.transAxes, fontsize=14, fontweight="bold")
    ax.text(0.5, 0.45,
            f"Found {n_rds} rlog RDS file(s) in deseq2/.\n"
            "Export VST counts to TSV to enable this plot.",
            ha="center", va="center", transform=ax.transAxes, fontsize=9, color="#555555")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print("PCA plot: placeholder written (no VST matrix available)")
