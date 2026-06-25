"""Merge CheckM2 quality + GTDB-Tk taxonomy for all samples into one TSV."""
import pandas as pd
import os

results_dir = snakemake.params.results_dir
samples     = snakemake.params.samples
out_path    = snakemake.output[0]

rows = []
for s in samples:
    qc_f = f"{results_dir}/{s}/checkm2/quality_report.tsv"
    tx_f = f"{results_dir}/{s}/gtdbtk/gtdbtk.bac120.summary.tsv"
    if not os.path.exists(qc_f):
        continue
    qc = pd.read_csv(qc_f, sep="\t")
    qc.insert(0, "sample", s)

    if os.path.exists(tx_f) and os.path.getsize(tx_f) > 0:
        try:
            tx = pd.read_csv(tx_f, sep="\t")[["user_genome", "classification"]]
            tx = tx.rename(columns={"user_genome": "Name", "classification": "gtdbtk_taxonomy"})
            qc = qc.merge(tx, on="Name", how="left")
        except Exception:
            pass

    rows.append(qc)

if not rows:
    pd.DataFrame().to_csv(out_path, sep="\t", index=False)
else:
    df = pd.concat(rows, ignore_index=True)
    df["quality_score"] = df["Completeness"] - 5 * df["Contamination"]
    df["quality_tier"]  = pd.cut(
        df["quality_score"],
        bins=[-999, 50, 90, 999],
        labels=["Medium-quality", "High-quality", "Near-complete"],
    )
    df.to_csv(out_path, sep="\t", index=False)
    print(f"MAG summary: {len(df)} MAGs from {len(rows)} samples")
