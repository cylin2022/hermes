"""PCA plot from VST count matrix."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

df  = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
out = snakemake.output[0]

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
