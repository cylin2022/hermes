"""Parse wgs_snp variant_stats.txt → structured TSV."""
import pandas as pd
from pathlib import Path

stats_path  = snakemake.input.stats
results_dir = snakemake.params.results_dir
out_path    = snakemake.output[0]

rows = []
current = {}
for line in Path(stats_path).read_text().splitlines():
    line = line.strip()
    if not line:
        continue
    if "\t" in line:
        k, v = line.split("\t", 1)
        current[k.strip()] = v.strip()

df = pd.DataFrame([current]) if current else pd.DataFrame()
df.to_csv(out_path, sep="\t", index=False)
print(f"Variant summary written: {len(df)} rows")
