# Snakemake script: Manhattan plots (GWAS + Fst), QQ plot, PCA plot
suppressPackageStartupMessages({
    library(data.table)
    library(CMplot)
    library(ggplot2)
    library(dplyr)
})

gwas_file   <- snakemake@input[["gwas_salt"]]
fst_snp_file<- snakemake@input[["fst_snp"]]
fst_win_file<- snakemake@input[["fst_win"]]
pca_file    <- snakemake@input[["pca_vec"]]
pheno_file  <- snakemake@input[["pheno"]]

mht_gwas    <- snakemake@output[["mht_gwas"]]
mht_fst     <- snakemake@output[["mht_fst"]]
qq_out      <- snakemake@output[["qq"]]
pca_out     <- snakemake@output[["pca_plot"]]

gwas_p      <- as.numeric(snakemake@params[["gwas_p"]])
suggestive  <- as.numeric(snakemake@params[["suggestive"]])
fst_top_pct <- as.numeric(snakemake@params[["fst_top"]])
outdir      <- snakemake@params[["outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── 0. Chromosome name mapping (NCBI accession → integer) ────────────────────
# Non-model organisms use NCBI accession IDs (NC_031965.2, etc.).
# CMplot requires integer chromosome numbers; we sort NC_ accessions
# lexicographically and assign 1..N. Unplaced scaffolds (NW_, MT) are dropped.
make_chr_map <- function(chr_vec) {
    nc <- sort(unique(grep("^NC_", chr_vec, value = TRUE)))
    if (length(nc) == 0) return(NULL)
    map <- data.table(orig = nc, Chr = seq_along(nc))
    map
}

# ── 1. GWAS Manhattan + QQ ────────────────────────────────────────────────────
message("[plots] Reading GEMMA output: ", gwas_file)
gwas <- fread(gwas_file)

# GEMMA .assoc.txt columns: chr rs ps n_miss allele1 allele0 af beta se logl_H1 l_remle p_wald p_lrt p_score
# CMplot expects: SNP, Chr, Pos, P
# Use p_lrt as primary (most powerful of the three GEMMA tests)
chr_map <- make_chr_map(gwas$chr)
if (!is.null(chr_map)) {
    gwas <- merge(gwas, chr_map, by.x = "chr", by.y = "orig", all.x = FALSE)
    gwas_plt <- gwas[, .(SNP = rs, Chr = Chr, Pos = ps, P = p_lrt)]
} else {
    gwas_plt <- gwas[, .(SNP = rs, Chr = chr, Pos = ps, P = p_lrt)]
}
gwas_plt <- gwas_plt[!is.na(P) & P > 0]

# Genomic inflation factor λ
chi2   <- qchisq(1 - gwas_plt$P, df = 1)
lambda <- round(median(chi2, na.rm = TRUE) / qchisq(0.5, 1), 3)
message(sprintf("[plots] Genomic inflation factor λ = %.3f", lambda))
writeLines(sprintf("lambda = %.3f", lambda), file.path(outdir, "lambda.txt"))

# Manhattan plot — CMplot always writes to cwd; setwd to outdir first
message("[plots] Drawing GWAS Manhattan plot")
old_wd <- getwd()
setwd(outdir)
CMplot(
    gwas_plt,
    type            = "p",
    plot.type       = "m",
    threshold       = c(gwas_p, suggestive),
    threshold.lwd   = c(1, 1),
    threshold.lty   = c(1, 2),
    threshold.col   = c("red", "blue"),
    amplify         = TRUE,
    signal.cex      = 1.5,
    signal.pch      = 19,
    signal.col      = "red",
    file            = "pdf",
    file.name       = "gwas_salt",
    dpi             = 300,
    file.output     = TRUE,
    verbose         = FALSE
)
setwd(old_wd)
# CMplot names file Rect_Manhtn.gwas_salt.pdf — rename to declared output
cand <- list.files(outdir, pattern = "Rect_Manhtn\\.gwas_salt\\.pdf", full.names = TRUE)
if (length(cand)) file.rename(cand[1], mht_gwas)

# QQ plot — write to outdir then rename
message("[plots] Drawing QQ plot")
setwd(outdir)
CMplot(
    gwas_plt,
    type         = "p",
    plot.type    = "q",
    conf.int     = TRUE,
    conf.int.col = "#00000033",
    file         = "pdf",
    file.name    = "qq_gwas",
    file.output  = TRUE,
    verbose      = FALSE
)
setwd(old_wd)
cand_qq <- list.files(outdir, pattern = "QQplot\\.qq_gwas\\.pdf", full.names = TRUE)
if (length(cand_qq)) file.rename(cand_qq[1], qq_out)

# ── 2. Fst Manhattan ─────────────────────────────────────────────────────────
message("[plots] Reading per-SNP Fst: ", fst_snp_file)
fst <- fread(fst_snp_file)
# VCFtools output: CHROM POS WEIR_AND_COCKERHAM_FST
setnames(fst, c("CHROM", "POS", "FST"))
fst <- fst[!is.na(FST) & FST >= 0]

fst_threshold <- quantile(fst$FST, fst_top_pct / 100, na.rm = TRUE)
message(sprintf("[plots] Fst top %g%% threshold = %.4f", 100 - fst_top_pct, fst_threshold))
writeLines(sprintf("fst_outlier_threshold = %.4f", fst_threshold), file.path(outdir, "fst_threshold.txt"))

fst_chr_map <- make_chr_map(fst$CHROM)
if (!is.null(fst_chr_map)) {
    fst <- merge(fst, fst_chr_map, by.x = "CHROM", by.y = "orig", all.x = FALSE)
    fst_plt <- fst[, .(SNP = paste0(CHROM, ":", POS), Chr = Chr, Pos = POS, FST = FST)]
} else {
    fst_plt <- fst[, .(SNP = paste0(CHROM, ":", POS), Chr = CHROM, Pos = POS, FST = FST)]
}

setwd(outdir)
CMplot(
    fst_plt,
    type            = "p",
    plot.type       = "m",
    threshold       = fst_threshold,
    threshold.lwd   = 1,
    threshold.lty   = 1,
    threshold.col   = "red",
    amplify         = TRUE,
    signal.col      = "red",
    ylab            = expression(F[ST]),
    file            = "pdf",
    file.name       = "fst_persnp",
    file.output     = TRUE,
    verbose         = FALSE
)
setwd(old_wd)
cand_fst <- list.files(outdir, pattern = "Rect_Manhtn\\.fst_persnp\\.pdf", full.names = TRUE)
if (length(cand_fst)) file.rename(cand_fst[1], mht_fst)

# ── 3. PCA plot ───────────────────────────────────────────────────────────────
message("[plots] Drawing PCA plot")
pca  <- fread(pca_file)
setnames(pca, c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2))))

pheno <- fread(pheno_file, header = TRUE)
setnames(pheno, old = c("IID", "PHENOTYPE"), new = c("IID", "phenotype"), skip_absent = TRUE)
pca   <- merge(pca, pheno[, .(IID, phenotype)], by = "IID", all.x = TRUE)
pca[, group := ifelse(phenotype == 1, "Salt-tolerant", "Salt-intolerant")]

p <- ggplot(pca, aes(PC1, PC2, colour = group, label = IID)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_colour_manual(values = c("Salt-tolerant" = "#e41a1c", "Salt-intolerant" = "#377eb8")) +
    theme_bw(base_size = 13) +
    labs(title = "PCA — population structure", colour = "Group",
         x = "PC1", y = "PC2")
ggsave(pca_out, p, width = 7, height = 5)

message("[plots] Done. Output: ", outdir)
