"""Copy/convert the single merged UMAP PDF produced by scrnaseq scanpy_pipeline.

The scrnaseq workflow produces scanpy/umap_clusters.pdf (a single merged file),
not per-sample umap_leiden.png files. This script converts the first page of the
PDF to a PNG for embedding in the HTML report.
"""
import os
import shutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

umap_pdf = snakemake.input[0]
out_path = snakemake.output[0]

if not os.path.exists(umap_pdf) or os.path.getsize(umap_pdf) == 0:
    # Produce placeholder
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.text(0.5, 0.5, "UMAP not available", ha="center", va="center",
            transform=ax.transAxes, fontsize=12)
    ax.axis("off")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"UMAP plot: placeholder written (source not found: {umap_pdf})")
else:
    # Try PDF→PNG conversion via pypdf + matplotlib, fallback to placeholder
    try:
        from pdf2image import convert_from_path
        images = convert_from_path(umap_pdf, dpi=150, first_page=1, last_page=1)
        images[0].save(out_path)
        print(f"UMAP plot: converted from {umap_pdf}")
    except Exception as e:
        # pdf2image not available or conversion failed — produce informative placeholder
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(0.5, 0.6, "UMAP clusters", ha="center", va="center",
                transform=ax.transAxes, fontsize=14, fontweight="bold")
        ax.text(0.5, 0.45,
                f"Source: {os.path.basename(umap_pdf)}\n"
                "(PDF→PNG conversion unavailable; install pdf2image + poppler)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9,
                color="#555555")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"UMAP plot: placeholder written (pdf2image error: {e})")
