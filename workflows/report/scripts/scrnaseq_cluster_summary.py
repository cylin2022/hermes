"""Read the single merged marker_genes.csv produced by the scrnaseq scanpy_pipeline rule."""
import pandas as pd
import os

results_dir = snakemake.params.results_dir
out_path    = snakemake.output[0]

# scrnaseq produces a single merged file: scanpy/marker_genes.csv
marker_file = snakemake.input[0]

if os.path.exists(marker_file) and os.path.getsize(marker_file) > 0:
    df = pd.read_csv(marker_file)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Cluster summary: {len(df)} marker rows from {marker_file}")
else:
    pd.DataFrame().to_csv(out_path, sep="\t", index=False)
    print(f"Cluster summary: empty (file not found or empty: {marker_file})")
