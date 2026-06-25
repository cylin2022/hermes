#!/usr/bin/env Rscript
# plot_results.R — Visualization for genomic prediction CV results
#
# Usage:
#   Rscript plot_results.R <cv_summary.csv> <cv_metrics_full.rds> <pca_eigenvec> <outdir>

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(pROC)
  library(data.table)
  library(dplyr)
  library(scales)
})

args        <- commandArgs(trailingOnly = TRUE)
summary_csv <- args[1]
metrics_rds <- args[2]
pca_file    <- args[3]   # PLINK2 .eigenvec
pheno_file  <- args[4]
outdir      <- args[5]

dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

summary_df   <- read.csv(summary_csv)
fold_results <- readRDS(metrics_rds)

MODEL_LABELS <- c(gblup = "GBLUP", lasso = "LASSO",
                  rf = "Random Forest", xgb = "XGBoost")
MODEL_COLORS <- c(GBLUP = "#4E79A7", LASSO = "#F28E2B",
                  "Random Forest" = "#59A14F", XGBoost = "#E15759")

# ── Plot 1: AUC bar chart with SD error bars ──────────────────────────────────
summary_df$model <- factor(summary_df$model,
                            levels = c("GBLUP", "LASSO", "Random Forest", "XGBoost"))
p_auc <- ggplot(summary_df, aes(x = model, y = mean_auc, fill = model)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
                width = 0.25, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", mean_auc)),
            vjust = -0.8, size = 4) +
  scale_fill_manual(values = MODEL_COLORS) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = percent_format(accuracy = 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  labs(title = "5-fold CV AUC (mean ± SD)",
       subtitle = "Salt tolerance genomic prediction — tilapia WGS",
       x = NULL, y = "AUC-ROC") +
  theme_cowplot(12) +
  theme(plot.subtitle = element_text(colour = "grey40"))

ggsave(file.path(outdir, "plots", "auc_comparison.pdf"),
       p_auc, width = 7, height = 5)

# ── Plot 2: All-metrics heatmap (AUC, Accuracy, Sensitivity, Specificity) ────
metric_long <- summary_df %>%
  select(model, mean_auc, mean_accuracy, mean_sens, mean_spec) %>%
  tidyr::pivot_longer(-model, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    mean_auc      = "AUC-ROC",
    mean_accuracy = "Accuracy",
    mean_sens     = "Sensitivity",
    mean_spec     = "Specificity"
  ))

p_heat <- ggplot(metric_long, aes(x = metric, y = model, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 4) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                        midpoint = 0.7, limits = c(0.4, 1.0),
                        name = "Score") +
  labs(title = "Performance metrics (5-fold CV mean)",
       x = NULL, y = NULL) +
  theme_cowplot(12) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank())

ggsave(file.path(outdir, "plots", "metrics_heatmap.pdf"),
       p_heat, width = 7, height = 4)

# ── Plot 3: ROC curves (aggregated predictions across all folds) ──────────────
roc_plots <- list()
models <- c("gblup", "lasso", "rf", "xgb")

roc_df <- lapply(models, function(m) {
  probs  <- unlist(lapply(fold_results, function(f) f[[m]]$prob))
  truths <- unlist(lapply(fold_results, function(f) f$truth))
  if (all(is.na(probs))) return(NULL)
  roc_obj <- roc(truths, probs, quiet = TRUE)
  data.frame(
    fpr   = 1 - roc_obj$specificities,
    tpr   = roc_obj$sensitivities,
    model = MODEL_LABELS[[m]],
    auc   = as.numeric(auc(roc_obj))
  )
})
roc_df <- do.call(rbind, Filter(Negate(is.null), roc_df))
roc_df$model_label <- sprintf("%s (AUC=%.3f)", roc_df$model, roc_df$auc)

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = model)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = MODEL_COLORS, name = NULL,
                      labels = unique(roc_df$model_label[order(roc_df$model)])) +
  labs(title = "ROC Curves (aggregated 5-fold CV)",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_cowplot(12) +
  theme(legend.position = c(0.55, 0.2))

ggsave(file.path(outdir, "plots", "roc_curves.pdf"),
       p_roc, width = 6, height = 6)

# ── Plot 4: Feature importance (LASSO + RF + XGBoost) ────────────────────────
TOP_N <- 20

# LASSO: average |coefficient| across folds (only SNPs selected in ≥1 fold)
lasso_coefs <- lapply(fold_results, function(f) {
  co <- f$lasso$coef
  if (is.null(co)) return(NULL)
  df <- data.frame(snp = rownames(co)[-1], coef = as.numeric(co[-1]))
  df[df$coef != 0, ]
})
lasso_coefs <- do.call(rbind, Filter(Negate(is.null), lasso_coefs))
if (!is.null(lasso_coefs) && nrow(lasso_coefs) > 0) {
  lasso_imp <- lasso_coefs %>%
    group_by(snp) %>%
    summarise(mean_abs_coef = mean(abs(coef)), n_folds_selected = n()) %>%
    arrange(desc(mean_abs_coef)) %>%
    head(TOP_N)

  p_lasso_imp <- ggplot(lasso_imp, aes(x = reorder(snp, mean_abs_coef),
                                        y = mean_abs_coef,
                                        fill = n_folds_selected)) +
    geom_col() +
    coord_flip() +
    scale_fill_viridis_c(name = "# folds\nselected", limits = c(1, 5)) +
    labs(title = sprintf("LASSO: Top %d SNPs by |coefficient|", TOP_N),
         subtitle = "Mean |β| across folds in which SNP was selected",
         x = NULL, y = "Mean |β|") +
    theme_cowplot(10)
  ggsave(file.path(outdir, "plots", "lasso_importance.pdf"),
         p_lasso_imp, width = 8, height = 6)
}

# RF: average impurity importance across folds
rf_imps <- lapply(fold_results, function(f) {
  imp <- f$rf$importance
  if (is.null(imp) || all(is.na(imp))) return(NULL)
  data.frame(snp = names(imp), importance = as.numeric(imp))
})
rf_imps <- do.call(rbind, Filter(Negate(is.null), rf_imps))
if (!is.null(rf_imps) && nrow(rf_imps) > 0) {
  rf_imp <- rf_imps %>%
    group_by(snp) %>%
    summarise(mean_importance = mean(importance)) %>%
    arrange(desc(mean_importance)) %>%
    head(TOP_N)

  p_rf_imp <- ggplot(rf_imp, aes(x = reorder(snp, mean_importance),
                                   y = mean_importance)) +
    geom_col(fill = MODEL_COLORS["Random Forest"]) +
    coord_flip() +
    labs(title = sprintf("Random Forest: Top %d SNPs (Gini impurity)", TOP_N),
         subtitle = "Mean importance across 5 folds",
         x = NULL, y = "Mean Gini Importance") +
    theme_cowplot(10)
  ggsave(file.path(outdir, "plots", "rf_importance.pdf"),
         p_rf_imp, width = 8, height = 6)
}

# XGBoost: average gain importance across folds
xgb_imps <- lapply(fold_results, function(f) {
  imp <- f$xgb$importance
  if (is.null(imp) || !is.data.frame(imp)) return(NULL)
  imp[, c("Feature", "Gain")]
})
xgb_imps <- do.call(rbind, Filter(Negate(is.null), xgb_imps))
if (!is.null(xgb_imps) && nrow(xgb_imps) > 0) {
  xgb_imp <- xgb_imps %>%
    group_by(Feature) %>%
    summarise(mean_gain = mean(Gain)) %>%
    arrange(desc(mean_gain)) %>%
    head(TOP_N)

  p_xgb_imp <- ggplot(xgb_imp, aes(x = reorder(Feature, mean_gain),
                                     y = mean_gain)) +
    geom_col(fill = MODEL_COLORS["XGBoost"]) +
    coord_flip() +
    labs(title = sprintf("XGBoost: Top %d SNPs (Gain)", TOP_N),
         subtitle = "Mean gain across 5 folds",
         x = NULL, y = "Mean Gain") +
    theme_cowplot(10)
  ggsave(file.path(outdir, "plots", "xgb_importance.pdf"),
         p_xgb_imp, width = 8, height = 6)
}

# ── Plot 5: PCA coloured by phenotype ─────────────────────────────────────────
if (file.exists(pca_file)) {
  pca_data <- fread(pca_file, data.table = FALSE)
  pheno    <- fread(pheno_file, header = FALSE, data.table = FALSE,
                    col.names = c("sample", "phenotype"))
  colnames(pca_data)[1:2] <- c("FID", "IID")
  pca_data <- merge(pca_data, pheno, by.x = "IID", by.y = "sample")
  pca_data$Group <- ifelse(pca_data$phenotype == 1,
                            "Salt-tolerant", "Salt-intolerant")

  pc_cols <- grep("^PC", colnames(pca_data), value = TRUE)
  p_pca <- ggplot(pca_data, aes_string(x = pc_cols[1], y = pc_cols[2],
                                        colour = "Group")) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linewidth = 0.7) +
    scale_colour_manual(values = c("Salt-tolerant"   = "#E15759",
                                   "Salt-intolerant" = "#4E79A7"),
                        name = NULL) +
    labs(title = "PCA of SNP genotypes",
         subtitle = "Coloured by salt tolerance phenotype",
         x = "PC1", y = "PC2") +
    theme_cowplot(12) +
    theme(legend.position = "top")

  ggsave(file.path(outdir, "plots", "pca_phenotype.pdf"),
         p_pca, width = 6, height = 6)
}

message("All plots written to: ", file.path(outdir, "plots"))
