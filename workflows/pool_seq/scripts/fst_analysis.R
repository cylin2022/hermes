library(data.table)
library(ggplot2)

# ── parameters ────────────────────────────────────────────────────────────────
samples    <- snakemake@params[["samples"]]
pool_sizes <- setNames(as.integer(snakemake@params[["pool_sizes"]]), samples)
ploidy     <- as.integer(snakemake@params[["ploidy"]])   # 2 for diploid
comparison <- snakemake@params[["comparison"]]
win_size   <- as.integer(snakemake@params[["window"]])
step_size  <- as.integer(snakemake@params[["step"]])
min_cov    <- snakemake@params[["min_cov"]]
max_cov    <- snakemake@params[["max_cov"]]

pools_a   <- comparison$pools_a
pools_b   <- comparison$pools_b
comp_name <- comparison$name

# Validate pool names against samplesheet
missing <- setdiff(c(pools_a, pools_b), samples)
if (length(missing) > 0)
  stop(sprintf("[%s] Unknown pool names: %s. Available: %s",
               comp_name, paste(missing, collapse = ", "), paste(samples, collapse = ", ")))

# Haploid count = number of fish × ploidy (Hudson Fst requires haploid allele count)
n_a <- sum(pool_sizes[pools_a]) * ploidy
n_b <- sum(pool_sizes[pools_b]) * ploidy
cat(sprintf("[%s] Group A: %s (n=%d haploid)  Group B: %s (n=%d haploid)\n",
            comp_name, paste(pools_a, collapse="+"), n_a,
            paste(pools_b, collapse="+"), n_b))

# ── read AD table ─────────────────────────────────────────────────────────────
# columns: chr pos ref alt  sample1_AD  sample2_AD ...
# AD field format: "ref_count,alt_count"
col_names <- c("chr", "pos", "ref", "alt", paste0(samples, "_AD"))
dt <- fread(snakemake@input[["ad"]], header = FALSE, col.names = col_names, sep = "\t")
cat(sprintf("[%s] SNPs in VCF: %d\n", comp_name, nrow(dt)))

# ── parse AD → frequency + depth per pool ────────────────────────────────────
for (s in samples) {
  ad_col <- paste0(s, "_AD")
  parts  <- strsplit(dt[[ad_col]], ",", fixed = TRUE)
  ref_ct <- as.integer(vapply(parts, `[[`, character(1), 1))
  alt_ct <- vapply(parts, function(x) if (length(x) >= 2) as.integer(x[2]) else 0L, integer(1))
  depth  <- ref_ct + alt_ct
  dt[, (paste0(s, "_freq"))  := ifelse(depth > 0, alt_ct / depth, NA_real_)]
  dt[, (paste0(s, "_depth")) := depth]
}

# ── coverage filter ───────────────────────────────────────────────────────────
depth_cols <- paste0(samples, "_depth")
keep <- rowSums(dt[, ..depth_cols] >= min_cov) == length(samples) &
        rowSums(dt[, ..depth_cols] <= max_cov) == length(samples)
dt   <- dt[keep]
cat(sprintf("[%s] SNPs after coverage filter: %d\n", comp_name, nrow(dt)))

if (nrow(dt) == 0)
  stop(sprintf("[%s] No SNPs passed coverage filter (min_cov=%d, max_cov=%d). Check BAM coverage.",
               comp_name, min_cov, max_cov))

# ── group-level allele frequencies (weighted mean by pool size) ───────────────
freq_a_cols <- paste0(pools_a, "_freq")
freq_b_cols <- paste0(pools_b, "_freq")
w_a <- pool_sizes[pools_a] / sum(pool_sizes[pools_a])
w_b <- pool_sizes[pools_b] / sum(pool_sizes[pools_b])

dt[, p_a := as.vector(as.matrix(.SD) %*% w_a), .SDcols = freq_a_cols]
dt[, p_b := as.vector(as.matrix(.SD) %*% w_b), .SDcols = freq_b_cols]

# ── Hudson Fst (Bhatia et al. 2013) ──────────────────────────────────────────
# Appropriate for Pool-seq: corrects for finite haploid pool size.
# n_a / n_b = haploid chromosome count (pool_size * ploidy), NOT fish count.
# Fst = [(p_a - p_b)^2 - p_a(1-p_a)/(n_a-1) - p_b(1-p_b)/(n_b-1)] /
#       [p_a(1-p_b) + p_b(1-p_a)]
dt[, fst := {
  num <- (p_a - p_b)^2 -
         p_a * (1 - p_a) / (n_a - 1) -
         p_b * (1 - p_b) / (n_b - 1)
  den <- p_a * (1 - p_b) + p_b * (1 - p_a)
  pmax(0, ifelse(den > 0, num / den, NA_real_))
}]

cat(sprintf("[%s] Mean per-SNP Fst: %.4f\n", comp_name,
            mean(dt$fst, na.rm = TRUE)))

# ── top 1% SNPs ───────────────────────────────────────────────────────────────
fst_99 <- quantile(dt$fst, 0.99, na.rm = TRUE)
top    <- dt[fst >= fst_99, .(chr, pos, ref, alt, p_a, p_b, fst)]
setorder(top, -fst)
fwrite(top, snakemake@output[["top_snps"]], sep = "\t")
cat(sprintf("[%s] Top 1%% SNPs: %d (Fst >= %.4f)\n", comp_name, nrow(top), fst_99))

# ── sliding window mean Fst ───────────────────────────────────────────────────
chrs <- unique(dt$chr)
windows <- rbindlist(lapply(chrs, function(chr_id) {
  sub <- dt[chr == chr_id & !is.na(fst)]
  if (nrow(sub) < 5) return(NULL)
  max_pos <- max(sub$pos)
  starts  <- seq(1L, max_pos, by = step_size)
  rbindlist(lapply(starts, function(s) {
    e <- s + win_size - 1L
    w <- sub[pos >= s & pos <= e]
    if (nrow(w) < 3) return(NULL)
    data.table(chr = chr_id, start = s, end = e,
               mid = (s + e) %/% 2L,
               fst_mean = mean(w$fst),
               n_snps   = nrow(w))
  }))
}))

if (is.null(windows) || nrow(windows) == 0)
  stop(sprintf("[%s] No sliding windows with >= 3 SNPs. Increase window_size or check coverage.", comp_name))

fwrite(windows, snakemake@output[["windows"]], sep = "\t")

# ── Manhattan plot ────────────────────────────────────────────────────────────
chr_order <- chrs
chr_lens  <- windows[, .(len = max(end)), by = chr]
setorder(chr_lens, chr)
offsets   <- setNames(
  cumsum(c(0L, chr_lens$len[-nrow(chr_lens)])),
  chr_lens$chr
)

windows[, x       := offsets[chr] + mid]
windows[, chr_idx := match(chr, chr_order)]

fst_win_99  <- quantile(windows$fst_mean, 0.99, na.rm = TRUE)
tick_pos    <- offsets[chr_lens$chr] + chr_lens$len / 2
# Shorten NC_ accessions to numeric part for readability
tick_labels <- sub("^NC_0+([0-9]+)\\.[0-9]+$", "\\1", chr_lens$chr)

p <- ggplot(windows, aes(x, fst_mean, color = factor(chr_idx %% 2))) +
  geom_point(size = 0.4, alpha = 0.7) +
  geom_hline(yintercept = fst_win_99,
             linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  scale_color_manual(
    values = c("0" = "#2166AC", "1" = "#92C5DE"),
    guide  = "none"
  ) +
  scale_x_continuous(
    breaks = tick_pos,
    labels = tick_labels,
    expand = c(0.01, 0)
  ) +
  labs(
    title    = sprintf("Pool-seq Fst — %s", comp_name),
    subtitle = sprintf(
      "Group A: %s  |  Group B: %s  |  Window: %g kb / step: %g kb  |  top 1%% threshold: %.3f",
      paste(pools_a, collapse = " + "),
      paste(pools_b, collapse = " + "),
      win_size  / 1000,
      step_size / 1000,
      fst_win_99
    ),
    x = "Chromosome",
    y = "Mean Fst (sliding window)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(snakemake@output[["pdf"]], p,
       width = 14, height = 5, device = grDevices::pdf)

cat(sprintf("[%s] Done.\n", comp_name))
