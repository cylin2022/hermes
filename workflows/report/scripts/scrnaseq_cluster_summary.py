"""Concatenate per-sample cluster_markers.tsv files across all samples."""
import pandas as pd
import os

results_dir = snakemake.params.results_dir
samples     = snakemake.params.samples
out_path    = snakemake.output[0]

frames = []
for s in samples:
    f = f"{results_dir}/{s}/scanpy/cluster_markers.tsv"
    if os.path.exists(f):
        df = pd.read_csv(f, sep="\t")
        df.insert(0, "sample", s)
        frames.append(df)

if frames:
    pd.concat(frames, ignore_index=True).to_csv(out_path, sep="\t", index=False)
else:
    pd.DataFrame().to_csv(out_path, sep="\t", index=False)
