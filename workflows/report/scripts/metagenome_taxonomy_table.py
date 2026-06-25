"""Merge Bracken species-level abundance tables across all samples."""
import pandas as pd
import os

results_dir = snakemake.params.results_dir
samples     = snakemake.params.samples
out_path    = snakemake.output[0]

frames = []
for s in samples:
    f = f"{results_dir}/{s}/bracken/{s}.S.bracken"
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t")
    df = df[["name", "fraction_total_reads"]].rename(
        columns={"name": "species", "fraction_total_reads": s}
    )
    frames.append(df.set_index("species"))

if not frames:
    pd.DataFrame().to_csv(out_path, sep="\t")
else:
    merged = pd.concat(frames, axis=1).fillna(0)
    merged["mean_abundance"] = merged.mean(axis=1)
    merged = merged.sort_values("mean_abundance", ascending=False).drop(columns="mean_abundance")
    merged.to_csv(out_path, sep="\t")
    print(f"Taxonomy table: {len(merged)} species, {len(frames)} samples")
