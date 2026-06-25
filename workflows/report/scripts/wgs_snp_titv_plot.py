"""Ti/Tv bar plot from variant summary TSV."""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df  = pd.read_csv(snakemake.input[0], sep="\t")
out = snakemake.output[0]

fig, ax = plt.subplots(figsize=(max(4, len(df) * 0.6 + 2), 4))
if "TiTv" in df.columns and not df.empty:
    ax.bar(df.index if "sample" not in df.columns else df["sample"],
           df["TiTv"].astype(float), color="#4393c3")
    ax.set_ylabel("Ti/Tv ratio")
    ax.set_xlabel("Sample")
    ax.axhline(2.0, ls="--", color="gray", lw=0.8, label="Expected ~2.0")
    ax.legend()
else:
    ax.text(0.5, 0.5, "Ti/Tv data not available",
            ha="center", va="center", transform=ax.transAxes)
ax.set_title("Transition / Transversion Ratio")
plt.tight_layout()
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
