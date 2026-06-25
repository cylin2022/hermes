"""Summarise DESeq2 all_results.tsv: count significant genes per contrast."""
import pandas as pd
from pathlib import Path

in_path  = snakemake.input[0]
out_path = snakemake.output[0]

df = pd.read_csv(in_path, sep="\t")
if "padj" in df.columns and "log2FoldChange" in df.columns:
    df["significant"] = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)
df.to_csv(out_path, sep="\t", index=False)
print(f"DE summary: {len(df)} genes")
