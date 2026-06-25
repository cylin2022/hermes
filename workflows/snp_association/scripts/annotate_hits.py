"""
Snakemake script: merge GWAS (salt adjusted + body size) + Fst → classify SNPs

SNP categories:
  salt_specific : p_salt < suggestive  AND  p_size >= suggestive (or no size GWAS)
  size_specific : p_size < suggestive  AND  p_salt >= suggestive
  pleiotropic   : both p_salt AND p_size < suggestive  → affects both traits
  fst_only      : Fst outlier but not GWAS-significant (population structure signal)

Finds nearest gene (±flank_bp) for each candidate using BEDTools.
"""

import subprocess
import os
from pathlib import Path
import pandas as pd
import numpy as np

gwas_salt    = snakemake.input["gwas_salt"]
gwas_size    = snakemake.input.get("gwas_size", None)
fst_file     = snakemake.input["fst_snp"]
annot_file   = snakemake.input["annot"]
cand_snps    = snakemake.output["cand_snps"]
cand_genes   = snakemake.output["cand_genes"]

gwas_p       = float(snakemake.params["gwas_p"])
suggestive   = float(snakemake.params["suggestive"])
fst_top_pct  = float(snakemake.params["fst_top_pct"])
flank_bp     = int(snakemake.params["flank_bp"])
has_cov      = bool(snakemake.params["has_cov"])
outdir       = snakemake.params["outdir"]

Path(outdir).mkdir(parents=True, exist_ok=True)

# ── Load GWAS results ─────────────────────────────────────────────────────────
def load_gemma(path, suffix):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    df = df.dropna(subset=["p_lrt"])
    df["chr_pos"] = df["chr"].astype(str) + ":" + df["ps"].astype(str)
    rename = {"p_lrt": f"p_lrt_{suffix}", "p_wald": f"p_wald_{suffix}",
              "beta": f"beta_{suffix}", "se": f"se_{suffix}", "af": f"af_{suffix}"}
    return df.rename(columns=rename)

salt = load_gemma(gwas_salt, "salt")

if has_cov and gwas_size and Path(gwas_size).exists() and Path(gwas_size).stat().st_size > 0:
    size = load_gemma(gwas_size, "size")
    merged = salt.merge(size[["chr_pos", "p_lrt_size", "p_wald_size", "beta_size"]],
                        on="chr_pos", how="left")
else:
    merged = salt.copy()
    merged["p_lrt_size"] = np.nan
    merged["beta_size"]  = np.nan
    has_cov = False

# ── Load Fst ─────────────────────────────────────────────────────────────────
fst = pd.read_csv(fst_file, sep="\t",
                  names=["CHROM", "POS", "FST"], skiprows=1)
fst = fst[fst["FST"] >= 0].dropna()
fst_threshold = np.percentile(fst["FST"], fst_top_pct)
fst["chr_pos"] = fst["CHROM"].astype(str) + ":" + fst["POS"].astype(str)
merged = merged.merge(fst[["chr_pos", "FST"]], on="chr_pos", how="left")

print(f"[annotate] Fst outlier threshold ({fst_top_pct}th pct): {fst_threshold:.4f}")

# ── Classify SNPs ─────────────────────────────────────────────────────────────
merged["sig_salt"] = merged["p_lrt_salt"] < suggestive
merged["sig_size"] = merged["p_lrt_size"] < suggestive if has_cov else False
merged["fst_out"]  = merged["FST"].fillna(0) >= fst_threshold

def classify(row):
    s, z, f = row["sig_salt"], row["sig_size"], row["fst_out"]
    if s and z:   return "pleiotropic"
    if s and f:   return "salt_fst"        # strongest: both GWAS and Fst
    if s:         return "salt_specific"
    if z and f:   return "size_fst"
    if z:         return "size_specific"
    if f:         return "fst_only"
    return None

merged["category"] = merged.apply(classify, axis=1)
candidates = merged[merged["category"].notna()].copy()
candidates = candidates.sort_values("p_lrt_salt")

print(f"\n[annotate] SNP classification summary:")
print(candidates["category"].value_counts().to_string())

# ── Prioritised candidates (exclude size-specific) ────────────────────────────
salt_candidates = candidates[~candidates["category"].isin(["size_specific", "size_fst"])]

# ── BEDTools: nearest gene ────────────────────────────────────────────────────
gene_bed = os.path.join(outdir, "genes.bed")
closest_results = {}

if annot_file and os.path.exists(annot_file):
    if annot_file.endswith(".gff") or annot_file.endswith(".gff3"):
        # GFF3: attributes use key=value; extract value after 'ID='
        awk = r"""awk '$3=="gene"{id=$9; gsub(/.*ID=/,"",id); gsub(/;.*/,"",id); printf "%s\t%d\t%d\t%s\t.\t%s\n",$1,$4-1,$5,id,$7}'"""
    else:
        # GTF: attributes use key "value"; pass full attribute field for parse_name
        awk = r"""awk '$3=="gene"{printf "%s\t%d\t%d\t%s\t.\t%s\n",$1,$4-1,$5,$9,$7}'"""
    os.system(f"grep -v '^#' {annot_file} | {awk} | sort -k1,1 -k2,2n > {gene_bed}")

    snp_bed = os.path.join(outdir, "candidates.bed")
    with open(snp_bed, "w") as f:
        for _, row in candidates.iterrows():
            pos = int(row["ps"])
            f.write(f"{row['chr']}\t{max(0, pos-flank_bp)}\t{pos+flank_bp}\t{row['rs']}\n")

    sorted_bed = snp_bed.replace(".bed", ".sorted.bed")
    os.system(f"sort -k1,1 -k2,2n {snp_bed} > {sorted_bed}")
    closest_out = os.path.join(outdir, "closest_genes.txt")
    subprocess.run(f"bedtools closest -a {sorted_bed} -b {gene_bed} -D ref > {closest_out}", shell=True)

    if os.path.exists(closest_out) and os.path.getsize(closest_out) > 0:
        cl = pd.read_csv(closest_out, sep="\t", header=None,
                         names=["qchr","qst","qen","snp_id","gchr","gst","gen",
                                "gene_info","score","strand","dist"])
        def parse_name(info):
            for tok in str(info).split(";"):
                if "gene_name" in tok: return tok.split('"')[1] if '"' in tok else tok.split("=")[-1].strip()
            for tok in str(info).split(";"):
                if "gene_id" in tok or "Name=" in tok:
                    return tok.split('"')[1] if '"' in tok else tok.split("=")[-1].strip()
            return "unknown"
        cl["gene_name"]   = cl["gene_info"].apply(parse_name)
        cl["distance_kb"] = (cl["dist"].abs() / 1000).round(1)
        candidates = candidates.merge(cl[["snp_id","gene_name","distance_kb"]],
                                      left_on="rs", right_on="snp_id", how="left").drop(columns="snp_id")
    else:
        candidates["gene_name"] = "N/A"
        candidates["distance_kb"] = "N/A"
else:
    candidates["gene_name"]   = "N/A"
    candidates["distance_kb"] = "N/A"

# ── Write outputs ─────────────────────────────────────────────────────────────
base_cols = ["chr","ps","rs","allele1","allele0","af_salt",
             "beta_salt","se_salt","p_lrt_salt","p_wald_salt",
             "beta_size","p_lrt_size",
             "FST","category","gene_name","distance_kb"]
out_cols = [c for c in base_cols if c in candidates.columns]
candidates[out_cols].to_csv(cand_snps, index=False)

# Gene summary — ranked by best salt-GWAS p-value, showing category breakdown
gene_summary = (
    candidates[candidates["gene_name"].notna() & (candidates["gene_name"] != "N/A")]
    .groupby("gene_name")
    .agg(
        best_p_salt    = ("p_lrt_salt",  "min"),
        best_fst       = ("FST",         "max"),
        n_snps         = ("rs",          "count"),
        categories     = ("category",    lambda x: "|".join(sorted(x.unique()))),
    )
    .reset_index()
    .sort_values("best_p_salt")
)
gene_summary.to_csv(cand_genes, index=False)

print(f"\n[annotate] Total candidate SNPs: {len(candidates)}")
print(f"[annotate] Salt-specific + salt_fst (true candidates): {len(candidates[candidates['category'].isin(['salt_specific','salt_fst'])])}")
print(f"[annotate] Pleiotropic (both traits): {len(candidates[candidates['category']=='pleiotropic'])}")
print(f"[annotate] Confounders removed (size-specific): {len(candidates[candidates['category'].isin(['size_specific','size_fst'])])}")
print(f"[annotate] Candidate genes: {len(gene_summary)}")
