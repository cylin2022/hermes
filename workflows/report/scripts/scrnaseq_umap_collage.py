"""Tile per-sample UMAP PNGs into a single collage image."""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math

results_dir = snakemake.params.results_dir
samples     = snakemake.params.samples
out_path    = snakemake.output[0]

img_paths = [
    f"{results_dir}/{s}/scanpy/umap_leiden.png"
    for s in samples
    if os.path.exists(f"{results_dir}/{s}/scanpy/umap_leiden.png")
]

if not img_paths:
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "No UMAP images found", ha="center", va="center", transform=ax.transAxes)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
else:
    n = len(img_paths)
    cols = min(n, 3)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4.5))
    axes = [axes] if n == 1 else axes.flat
    for ax, p, s in zip(axes, img_paths, samples):
        img = mpimg.imread(p)
        ax.imshow(img)
        ax.set_title(s, fontsize=9)
        ax.axis("off")
    for ax in list(axes)[n:]:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
