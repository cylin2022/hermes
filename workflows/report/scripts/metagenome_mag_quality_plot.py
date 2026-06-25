"""Scatter plot: MAG completeness vs contamination, coloured by quality tier."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(snakemake.input[0], sep="\t")
out = snakemake.output[0]

if df.empty:
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, "No MAGs found", ha="center", va="center", transform=plt.gca().transAxes)
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
else:
    palette = {"Near-complete": "#2166ac", "High-quality": "#74add1", "Medium-quality": "#fdae61"}
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=df, x="Contamination", y="Completeness",
        hue="quality_tier", palette=palette,
        s=60, alpha=0.8, ax=ax,
    )
    ax.axhline(90, ls="--", lw=0.8, color="gray")
    ax.axhline(50, ls=":",  lw=0.8, color="gray")
    ax.axvline(5,  ls="--", lw=0.8, color="gray")
    ax.axvline(10, ls=":",  lw=0.8, color="gray")
    ax.set_xlim(-0.5, max(df["Contamination"].max() * 1.1, 12))
    ax.set_ylim(-2, 102)
    ax.set_xlabel("Contamination (%)")
    ax.set_ylabel("Completeness (%)")
    ax.set_title(f"MAG Quality  (n={len(df)})")
    ax.legend(title="Tier", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
