"""Stacked bar chart of top-N species abundance across samples."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

df  = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
out = snakemake.output[0]
top_n = snakemake.params.top_n

if df.empty:
    plt.figure(figsize=(6, 4))
    plt.text(0.5, 0.5, "No taxonomy data", ha="center", va="center", transform=plt.gca().transAxes)
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
else:
    # Keep top N by mean abundance; lump rest as "Other"
    top_sp = df.mean(axis=1).nlargest(top_n).index.tolist()
    plot_df = df.loc[top_sp].T.copy()
    other = df.drop(index=top_sp).sum()
    plot_df["Other"] = other.values
    plot_df = plot_df * 100  # fraction → percent

    colors = sns.color_palette("tab20", len(plot_df.columns))
    fig, ax = plt.subplots(figsize=(max(8, len(plot_df) * 0.9), 6))
    plot_df.plot(kind="bar", stacked=True, ax=ax, color=colors, width=0.75)
    ax.set_ylabel("Relative Abundance (%)")
    ax.set_xlabel("Sample")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.legend(title="Species", bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=7, ncol=1)
    ax.set_title(f"Top-{top_n} Species Abundance")
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
