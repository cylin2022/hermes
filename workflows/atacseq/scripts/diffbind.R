# DiffBind differential chromatin accessibility analysis
# Called as a Snakemake script; outputs per-contrast results + volcano plot.
suppressPackageStartupMessages({
    library(DiffBind)
    library(ggplot2)
    library(BiocParallel)
})

register(MulticoreParam(snakemake@threads))

meta_file    <- snakemake@input[["metadata"]]
contrast_str <- snakemake@params[["contrast"]]
comparisons  <- snakemake@params[["comparisons"]]
padj_thresh  <- as.numeric(snakemake@params[["padj"]])
lfc_thresh   <- as.numeric(snakemake@params[["lfc"]])
outdir       <- snakemake@params[["outdir"]]
bam_dir      <- snakemake@params[["bam_dir"]]
peak_dir     <- snakemake@params[["peak_dir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

parts <- strsplit(contrast_str, "_vs_")[[1]]
cond1 <- parts[1]
cond2 <- parts[2]

meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Build DiffBind sample sheet
db_sheet <- data.frame(
    SampleID  = meta$sample,
    Condition = meta$group,
    bamReads  = file.path(bam_dir,  paste0(meta$sample, ".markdup.bam")),
    Peaks     = file.path(peak_dir, paste0(meta$sample, "_peaks.narrowPeak")),
    PeakCaller = "narrow",
    stringsAsFactors = FALSE
)

dba <- dba(sampleSheet = db_sheet)
dba <- dba.count(dba, bParallel = TRUE)

# PCA on all samples (written once, not per contrast)
pca_file <- file.path(outdir, "pca.pdf")
if (!file.exists(pca_file)) {
    pdf(pca_file, width = 6, height = 5)
    dba.plotPCA(dba, label = DBA_ID)
    dev.off()
}

# Contrast
dba <- dba.contrast(dba,
    group1 = dba$masks[[cond1]],
    group2 = dba$masks[[cond2]],
    name1  = cond1,
    name2  = cond2)
dba <- dba.analyze(dba, bParallel = TRUE)

res <- dba.report(dba, contrast = 1, th = 1, fold = 0)  # retrieve all, filter below
res_df <- as.data.frame(res)
colnames(res_df)[colnames(res_df) == "FDR"] <- "padj"
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.table(res_df, snakemake@output[["results"]], sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano plot
sig <- res_df[!is.na(res_df$padj) &
              res_df$padj < padj_thresh &
              abs(res_df$Fold) >= lfc_thresh, ]
p <- ggplot(res_df, aes(x = Fold, y = -log10(p.value))) +
    geom_point(alpha = 0.3, size = 0.7, colour = "grey60") +
    geom_point(data = sig, alpha = 0.8, size = 1.0, colour = "steelblue") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", colour = "firebrick") +
    geom_hline(yintercept = -log10(padj_thresh),         linetype = "dashed", colour = "firebrick") +
    labs(title = contrast_str, x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_bw(base_size = 12)
ggsave(snakemake@output[["volcano"]], p, width = 7, height = 5)
