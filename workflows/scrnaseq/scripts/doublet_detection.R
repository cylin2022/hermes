# Snakemake script: doublet detection with scDblFinder
suppressPackageStartupMessages({
    library(scDblFinder)
    library(SingleCellExperiment)
    library(BiocParallel)
})

mtx_dir  <- snakemake@input[["mtx_dir"]]
out_file <- snakemake@output[["scores"]]
sample_id <- snakemake@params[["sample"]]

message("[scDblFinder] Loading count matrix: ", mtx_dir)
sce <- read10xCounts(mtx_dir, col.names = TRUE)
colnames(sce) <- paste0(sample_id, "_", colnames(sce))

message("[scDblFinder] Running doublet detection on ", ncol(sce), " cells")
set.seed(42)
sce <- scDblFinder(sce, BPPARAM = SerialParam())

df <- data.frame(
    barcode        = colnames(sce),
    doublet_score  = sce$scDblFinder.score,
    doublet_class  = sce$scDblFinder.class,
    stringsAsFactors = FALSE
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write.csv(df, out_file, row.names = FALSE)
message("[scDblFinder] Doublets found: ", sum(df$doublet_class == "doublet"),
        " / ", nrow(df), " cells (", round(mean(df$doublet_class == "doublet") * 100, 1), "%)")
