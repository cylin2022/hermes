"""
Snakemake script: parse metadata CSV → PLINK pheno, GEMMA covariate,
VCFtools pop files, and body-size phenotype file.

Metadata CSV format (comma-separated, with header):
    sample_id, phenotype, [weight_g, length_cm, ...]
    fish_001,  1,          52.3,     16.1
    fish_002,  0,          41.7,     14.5

phenotype: 1 = salt-tolerant (case), 0 = salt-intolerant (control), NA = missing
Covariate columns (anything beyond the first two) are OPTIONAL.
When present:
  - Values are z-score standardised (mean=0, sd=1)
  - GEMMA covariate file is written (no header; first column = intercept=1)
  - Body-size phenotype file is written (first covariate column) for secondary GWAS
"""

import csv
import subprocess
import numpy as np
from pathlib import Path

metadata_file    = snakemake.input["metadata"]
vcf_file         = snakemake.input["vcf"]

plink_pheno_out  = snakemake.output["plink_pheno"]
pop_tol_out      = snakemake.output["pop_tolerant"]
pop_int_out      = snakemake.output["pop_intolerant"]
gemma_cov_out    = snakemake.output["gemma_cov"]
bodysize_pheno   = snakemake.output["bodysize_pheno"]

Path(plink_pheno_out).parent.mkdir(parents=True, exist_ok=True)

# ── Read metadata CSV ─────────────────────────────────────────────────────────
pheno    = {}
cov_data = {}
cov_cols = []

with open(metadata_file) as fh:
    reader = csv.DictReader(fh)
    cov_cols = [c for c in reader.fieldnames if c not in ("sample_id", "phenotype")]
    for row in reader:
        sid = row["sample_id"].strip()
        pval = row["phenotype"].strip()
        pheno[sid] = int(pval) if pval not in ("NA", "na", "", "N/A") else -9
        if cov_cols:
            cov_data[sid] = [
                float(row[c]) if row[c].strip() not in ("NA", "na", "", "N/A") else np.nan
                for c in cov_cols
            ]

has_cov = bool(cov_cols)

# ── Get sample order from VCF header ─────────────────────────────────────────
result = subprocess.run(
    ["bcftools", "query", "-l", vcf_file],
    capture_output=True, text=True, check=True
)
vcf_samples = [s.strip() for s in result.stdout.splitlines() if s.strip()]

missing = [s for s in pheno if s not in vcf_samples]
if missing:
    print(f"[prepare] WARNING: {len(missing)} metadata samples absent from VCF: {missing[:5]}")

# ── PLINK2 phenotype file (FID IID PHENOTYPE; --1 flag: 0=control, 1=case) ───
with open(plink_pheno_out, "w") as f:
    f.write("FID\tIID\tPHENOTYPE\n")
    for sid in vcf_samples:
        pval = pheno.get(sid, -9)
        f.write(f"{sid}\t{sid}\t{pval}\n")

# ── VCFtools population files ─────────────────────────────────────────────────
with open(pop_tol_out, "w") as ft, open(pop_int_out, "w") as fi:
    for sid in vcf_samples:
        if sid not in pheno or pheno[sid] == -9:
            continue
        (ft if pheno[sid] == 1 else fi).write(sid + "\n")

n_tol = sum(1 for v in pheno.values() if v == 1)
n_int = sum(1 for v in pheno.values() if v == 0)
print(f"[prepare] Salt-tolerant (case): {n_tol}  |  Intolerant (control): {n_int}")

# ── Covariates (optional) ─────────────────────────────────────────────────────
if has_cov:
    cov_matrix = []
    for sid in vcf_samples:
        if sid in cov_data:
            cov_matrix.append(cov_data[sid])
        else:
            cov_matrix.append([np.nan] * len(cov_cols))

    cov_arr = np.array(cov_matrix, dtype=float)
    means   = np.nanmean(cov_arr, axis=0)
    stds    = np.nanstd(cov_arr, axis=0)
    stds[stds == 0] = 1.0
    cov_z   = (cov_arr - means) / stds

    # GEMMA covariate file: no header; first column = intercept (1)
    Path(gemma_cov_out).parent.mkdir(parents=True, exist_ok=True)
    with open(gemma_cov_out, "w") as f:
        for row in cov_z:
            vals = ["1"] + ["NA" if np.isnan(v) else f"{v:.6f}" for v in row]
            f.write("\t".join(vals) + "\n")
    print(f"[prepare] GEMMA covariates written ({len(cov_cols)} columns: {cov_cols})")

    # Body-size phenotype file: first covariate column (z-scored weight)
    Path(bodysize_pheno).parent.mkdir(parents=True, exist_ok=True)
    with open(bodysize_pheno, "w") as f:
        for row in cov_z:
            val = row[0]
            f.write("NA\n" if np.isnan(val) else f"{val:.6f}\n")
    print(f"[prepare] Body-size phenotype written (column: {cov_cols[0]})")
else:
    # Write empty placeholder files so Snakemake output requirements are satisfied
    for out in [gemma_cov_out, bodysize_pheno]:
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).touch()
    print("[prepare] No covariate columns found; covariate adjustment will be skipped.")
