# DESeq2 differential expression + clusterProfiler KEGG enrichment
# Called as a Snakemake script; outputs per-contrast results, volcano, enrichment table.
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(clusterProfiler)
    library(BiocParallel)
    library(data.table)
})

register(MulticoreParam(snakemake@threads))

counts_file  <- snakemake@input[["counts"]]
meta_file    <- snakemake@input[["metadata"]]
contrast_str <- snakemake@params[["contrast"]]
padj_thresh  <- as.numeric(snakemake@params[["padj"]])
lfc_thresh   <- as.numeric(snakemake@params[["lfc"]])
species_name <- snakemake@params[["species"]]
outdir       <- snakemake@params[["outdir"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

parts <- strsplit(contrast_str, "_vs_")[[1]]
cond1 <- parts[1]   # numerator
cond2 <- parts[2]   # denominator / reference

# featureCounts output: first line is comment; columns 1–6 are annotation, 7+ are BAM paths
raw      <- fread(counts_file, skip = 1)
gene_ids <- raw[[1]]
bam_cols <- colnames(raw)[7:ncol(raw)]
cnt_mat  <- as.matrix(raw[, 7:ncol(raw), with = FALSE])
rownames(cnt_mat) <- gene_ids
# Recover sample names from BAM path: .../star/{sample}/Aligned...bam
colnames(cnt_mat) <- basename(dirname(bam_cols))

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample

keep    <- meta$sample[meta$group %in% c(cond1, cond2)]
meta_s  <- meta[keep, , drop = FALSE]
cnt_s   <- cnt_mat[, keep, drop = FALSE]
meta_s$group <- factor(meta_s$group, levels = c(cond2, cond1))

dds <- DESeqDataSetFromMatrix(countData = cnt_s, colData = meta_s, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds, parallel = TRUE)

res    <- results(dds, contrast = c("group", cond1, cond2), alpha = padj_thresh)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.table(res_df, snakemake@output[["results"]], sep = "\t", quote = FALSE, row.names = FALSE)

rld <- rlog(dds, blind = FALSE)
saveRDS(rld, snakemake@output[["rlog"]])

# Volcano plot
sig <- res_df[!is.na(res_df$padj) &
              res_df$padj < padj_thresh &
              abs(res_df$log2FoldChange) >= lfc_thresh, ]
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(alpha = 0.3, size = 0.7, colour = "grey60") +
    geom_point(data = sig, alpha = 0.8, size = 1.0, colour = "firebrick") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", colour = "steelblue") +
    geom_hline(yintercept = -log10(padj_thresh),         linetype = "dashed", colour = "steelblue") +
    labs(title = contrast_str, x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_bw(base_size = 12)
ggsave(snakemake@output[["volcano"]], p, width = 7, height = 5)

# KEGG enrichment — non-fatal if organism not mapped or too few genes
kegg_map <- c(
    Mus_musculus          = "mmu",
    Homo_sapiens          = "hsa",
    Danio_rerio           = "dre",
    Oreochromis_niloticus = "oni",
    Gallus_gallus         = "gga",
    Rattus_norvegicus     = "rno"
)
kegg_org  <- kegg_map[species_name]
sig_genes <- sig$gene[!is.na(sig$gene)]

enrich_df <- tryCatch({
    if (!is.na(kegg_org) && length(sig_genes) >= 5) {
        ego <- enrichKEGG(gene = sig_genes, organism = kegg_org, pvalueCutoff = 0.05)
        as.data.frame(ego)
    } else {
        data.frame()
    }
}, error = function(e) {
    message("clusterProfiler: ", conditionMessage(e))
    data.frame()
})
write.table(enrich_df, snakemake@output[["enrich"]], sep = "\t", quote = FALSE, row.names = FALSE)
