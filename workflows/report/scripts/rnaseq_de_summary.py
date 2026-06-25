"""Summarise DESeq2 per-contrast results: count significant genes per contrast.

The rnaseq workflow writes {contrast}_results.tsv files (not a single all_results.tsv).
This script globs all *_results.tsv files from the deseq2 directory and concatenates them,
adding a 'contrast' column derived from the filename.
"""
import pandas as pd
from pathlib import Path
import warnings

deseq2_dir = Path(snakemake.params.deseq2_dir)
out_path   = snakemake.output[0]

# Collect all per-contrast results files
result_files = sorted(deseq2_dir.glob("*_results.tsv"))

if not result_files:
    warnings.warn(f"No *_results.tsv files found in {deseq2_dir}. Writing empty summary.")
    pd.DataFrame().to_csv(out_path, sep="\t", index=False)
    print("DE summary: 0 contrasts found")
else:
    frames = []
    for f in result_files:
        contrast = f.stem.replace("_results", "")
        df = pd.read_csv(f, sep="\t")
        df.insert(0, "contrast", contrast)
        if "padj" in df.columns and "log2FoldChange" in df.columns:
            df["significant"] = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)
        frames.append(df)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(out_path, sep="\t", index=False)
    print(f"DE summary: {len(combined)} genes across {len(frames)} contrast(s)")
