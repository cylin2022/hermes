"""Parse wgs_snp stats/variant_summary.txt → structured TSV.

The file is produced by the wgs_snp variant_stats rule and contains sections:
  === PASS SNPs ===
  SN\t0\t<key>:\t<value>     (bcftools stats summary numbers)
  === Ts/Tv ratio (SNPs) ===
  TSTV\t0\t<ts>\t<tv>\t<ts/tv>\t<ts_1st_ALT>\t<tv_1st_ALT>\t<ts/tv_1st_ALT>
  === PASS INDELs ===
  SN\t0\t<key>:\t<value>
  === Per-sample SNP counts ===
  PSC\t0\t<id>\t<sample>\t<hom_RR>\t<hom_AA>\t<ts>\t<tv>\t<indels>...
"""
import pandas as pd
from pathlib import Path

stats_path  = snakemake.input.stats
results_dir = snakemake.params.results_dir
out_path    = snakemake.output[0]

rows = []
current = {}
section = ""

for line in Path(stats_path).read_text().splitlines():
    line = line.strip()
    if not line:
        continue

    # Track which section we are in
    if line.startswith("==="):
        section = line.strip("= ").strip()
        continue

    if not "\t" in line:
        continue

    parts = line.split("\t")
    tag = parts[0]

    if tag == "SN":
        # SN\t<id>\t<key>:\t<value>
        if len(parts) >= 4:
            key = parts[2].rstrip(":")
            val = parts[3]
            current[f"{section}|{key}"] = val

    elif tag == "TSTV":
        # TSTV\t<id>\t<ts>\t<tv>\t<ts/tv>\t...
        # Index 4 (0-based) is the Ts/Tv ratio
        if len(parts) >= 5:
            current[f"{section}|Ts/Tv ratio"] = parts[4]

if current:
    df = pd.DataFrame([current])
else:
    df = pd.DataFrame()

df.to_csv(out_path, sep="\t", index=False)
print(f"Variant summary written: {len(df)} rows, {len(df.columns)} columns")
