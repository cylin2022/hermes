#!/usr/bin/env Rscript
# run_cv_models.R
# 5-fold stratified CV genomic prediction: GBLUP, LASSO, Random Forest, XGBoost
#
# Usage:
#   Rscript run_cv_models.R <raw_file> <pheno_csv> <n_folds> <seed> <threads> <outdir>
#
# raw_file  : PLINK .raw (additive coded, 0/1/2, from plink2 --recode A)
# pheno_csv : sample_id,phenotype  (1=salt-tolerant, 0=salt-intolerant)
# outdir    : directory for cv_summary.csv and cv_metrics_full.rds

suppressPackageStartupMessages({
  library(rrBLUP)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(pROC)
  library(data.table)
})

args       <- commandArgs(trailingOnly = TRUE)
raw_file   <- args[1]
pheno_file <- args[2]
n_folds    <- as.integer(args[3])
seed       <- as.integer(args[4])
n_threads  <- as.integer(args[5])
outdir     <- args[6]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading SNP matrix...")
raw  <- fread(raw_file, data.table = FALSE)
sids <- raw$IID
# Columns 1-6: FID IID PAT MAT SEX PHENOTYPE; 7+ are SNP_allele columns
X    <- as.matrix(raw[, 7:ncol(raw)])
rownames(X) <- sids

# Impute missing values with column mean (mean imputation for missing genotypes)
for (j in seq_len(ncol(X))) {
  na_idx <- is.na(X[, j])
  if (any(na_idx)) X[na_idx, j] <- mean(X[!na_idx, j])
}

pheno <- fread(pheno_file, header = TRUE, data.table = FALSE)
colnames(pheno)[1:2] <- c("sample", "phenotype")
pheno <- pheno[match(sids, pheno$sample), ]
stopifnot(!anyNA(pheno$sample))      # all VCF samples must be in phenotype file
stopifnot(!anyNA(pheno$phenotype))   # phenotype must be non-missing for all samples
y     <- as.integer(pheno$phenotype)
n     <- length(y)

message(sprintf("Loaded %d samples × %d SNPs. Phenotype: %d cases, %d controls.",
                n, ncol(X), sum(y == 1), sum(y == 0)))

# ── Stratified 5-fold CV assignment ──────────────────────────────────────────
# Assign folds within each class so each fold has equal class balance
folds             <- integer(n)
folds[y == 1]     <- sample(rep(seq_len(n_folds), length.out = sum(y == 1)))
folds[y == 0]     <- sample(rep(seq_len(n_folds), length.out = sum(y == 0)))

# ── Metric helper ─────────────────────────────────────────────────────────────
binary_metrics <- function(prob, truth, threshold = 0.5) {
  pred <- as.integer(prob >= threshold)
  tp   <- sum(pred == 1 & truth == 1)
  tn   <- sum(pred == 0 & truth == 0)
  fp   <- sum(pred == 1 & truth == 0)
  fn   <- sum(pred == 0 & truth == 1)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  acc  <- (tp + tn) / length(truth)
  list(accuracy = acc, sensitivity = sens, specificity = spec)
}

# ── Per-fold storage ──────────────────────────────────────────────────────────
fold_results <- vector("list", n_folds)

for (fold in seq_len(n_folds)) {
  message(sprintf("\n=== Fold %d/%d ===", fold, n_folds))

  tr  <- which(folds != fold)
  te  <- which(folds == fold)
  Xtr <- X[tr, , drop = FALSE]
  Xte <- X[te, , drop = FALSE]
  ytr <- y[tr]
  yte <- y[te]

  # Scale using training-set statistics only (prevent data leakage)
  mu  <- colMeans(Xtr)
  sg  <- apply(Xtr, 2, sd)
  sg[sg == 0] <- 1
  Xtr_sc <- scale(Xtr, center = mu, scale = sg)
  Xte_sc <- scale(Xte, center = mu, scale = sg)

  fold_res <- list(fold = fold, truth = yte)

  # ── 1. GBLUP ────────────────────────────────────────────────────────────────
  # Use full dataset GRM with test phenotypes masked (standard genomic selection)
  message("  Running GBLUP...")
  tryCatch({
    X_all  <- rbind(Xtr_sc, Xte_sc)
    K_all  <- A.mat(X_all)
    ids    <- seq_len(nrow(X_all))

    df_gblup <- data.frame(
      id        = ids,
      phenotype = c(ytr, rep(NA, length(yte)))
    )
    fit_gblup <- kin.blup(data = df_gblup, geno = "id", pheno = "phenotype",
                          REML = TRUE, K = K_all)

    # fit_gblup$pred gives BLUP for all, including NAs (test individuals)
    te_ids       <- (length(tr) + 1):nrow(X_all)
    gebv_te      <- fit_gblup$pred[te_ids]
    prob_gblup   <- plogis(gebv_te)  # logistic link to [0,1]

    roc_obj      <- roc(yte, prob_gblup, quiet = TRUE)
    m            <- binary_metrics(prob_gblup, yte)
    fold_res$gblup <- c(auc = as.numeric(auc(roc_obj)), m,
                        list(prob = prob_gblup))
    message(sprintf("    AUC=%.3f  Acc=%.3f  Sens=%.3f  Spec=%.3f",
                    fold_res$gblup$auc, fold_res$gblup$accuracy,
                    fold_res$gblup$sensitivity, fold_res$gblup$specificity))
  }, error = function(e) {
    message("    GBLUP error: ", conditionMessage(e))
    fold_res$gblup <<- list(auc = NA, accuracy = NA, sensitivity = NA,
                             specificity = NA, prob = NA)
  })

  # ── 2. LASSO logistic regression ────────────────────────────────────────────
  message("  Running LASSO...")
  tryCatch({
    # Inner 5-fold CV on training set to select lambda
    cv_lasso <- cv.glmnet(Xtr_sc, ytr, family = "binomial",
                          alpha = 1, nfolds = 5, type.measure = "auc",
                          parallel = FALSE)
    prob_lasso <- as.numeric(
      predict(cv_lasso, Xte_sc, s = "lambda.1se", type = "response")
    )
    roc_obj    <- roc(yte, prob_lasso, quiet = TRUE)
    m          <- binary_metrics(prob_lasso, yte)
    coef_mat   <- coef(cv_lasso, s = "lambda.1se")
    n_selected <- sum(coef_mat[-1] != 0)  # exclude intercept

    fold_res$lasso <- c(auc = as.numeric(auc(roc_obj)), m,
                        n_selected = n_selected,
                        list(prob = prob_lasso,
                             coef = coef_mat))
    message(sprintf("    AUC=%.3f  Acc=%.3f  Sens=%.3f  Spec=%.3f  SNPs=%d",
                    fold_res$lasso$auc, fold_res$lasso$accuracy,
                    fold_res$lasso$sensitivity, fold_res$lasso$specificity,
                    n_selected))
  }, error = function(e) {
    message("    LASSO error: ", conditionMessage(e))
    fold_res$lasso <<- list(auc = NA, accuracy = NA, sensitivity = NA,
                             specificity = NA, n_selected = NA, prob = NA)
  })

  # ── 3. Random Forest ─────────────────────────────────────────────────────────
  message("  Running Random Forest...")
  tryCatch({
    df_tr <- data.frame(y = factor(ytr, levels = c(0, 1)), Xtr_sc)
    df_te <- data.frame(Xte_sc)

    fit_rf <- ranger(
      y ~ ., data = df_tr,
      num.trees       = 500,
      probability     = TRUE,
      importance      = "impurity",
      min.node.size   = 5,
      num.threads     = n_threads,
      seed            = seed + fold
    )
    prob_rf  <- predict(fit_rf, df_te)$predictions[, "1"]
    roc_obj  <- roc(yte, prob_rf, quiet = TRUE)
    m        <- binary_metrics(prob_rf, yte)

    fold_res$rf <- c(auc = as.numeric(auc(roc_obj)), m,
                     list(prob = prob_rf,
                          importance = sort(fit_rf$variable.importance,
                                           decreasing = TRUE)))
    message(sprintf("    AUC=%.3f  Acc=%.3f  Sens=%.3f  Spec=%.3f",
                    fold_res$rf$auc, fold_res$rf$accuracy,
                    fold_res$rf$sensitivity, fold_res$rf$specificity))
  }, error = function(e) {
    message("    RF error: ", conditionMessage(e))
    fold_res$rf <<- list(auc = NA, accuracy = NA, sensitivity = NA,
                          specificity = NA, prob = NA, importance = NA)
  })

  # ── 4. XGBoost ───────────────────────────────────────────────────────────────
  message("  Running XGBoost...")
  tryCatch({
    dtrain <- xgb.DMatrix(Xtr_sc, label = ytr)
    dtest  <- xgb.DMatrix(Xte_sc)

    params <- list(
      objective        = "binary:logistic",
      eta              = 0.05,
      max_depth        = 4,
      subsample        = 0.8,
      colsample_bytree = 0.8,
      eval_metric      = "auc",
      nthread          = n_threads
    )
    # Inner CV to find optimal nrounds
    cv_xgb <- xgb.cv(
      params = params, data = dtrain,
      nrounds = 500, nfold = 5,
      early_stopping_rounds = 30,
      verbose = FALSE
    )
    best_n <- cv_xgb$best_iteration

    fit_xgb  <- xgboost(params = params, data = dtrain,
                         nrounds = best_n, verbose = 0)
    prob_xgb <- predict(fit_xgb, dtest)
    roc_obj  <- roc(yte, prob_xgb, quiet = TRUE)
    m        <- binary_metrics(prob_xgb, yte)
    imp_xgb  <- xgb.importance(model = fit_xgb)

    fold_res$xgb <- c(auc = as.numeric(auc(roc_obj)), m,
                      nrounds = best_n,
                      list(prob = prob_xgb, importance = imp_xgb))
    message(sprintf("    AUC=%.3f  Acc=%.3f  Sens=%.3f  Spec=%.3f  rounds=%d",
                    fold_res$xgb$auc, fold_res$xgb$accuracy,
                    fold_res$xgb$sensitivity, fold_res$xgb$specificity, best_n))
  }, error = function(e) {
    message("    XGBoost error: ", conditionMessage(e))
    fold_res$xgb <<- list(auc = NA, accuracy = NA, sensitivity = NA,
                           specificity = NA, nrounds = NA, prob = NA,
                           importance = NA)
  })

  fold_results[[fold]] <- fold_res
}

# ── Aggregate summary ─────────────────────────────────────────────────────────
get_metric <- function(model_key, metric_key) {
  sapply(fold_results, function(f) {
    val <- f[[model_key]][[metric_key]]
    if (is.null(val)) NA else val
  })
}

models <- c("gblup", "lasso", "rf", "xgb")
labels <- c("GBLUP", "LASSO", "Random Forest", "XGBoost")

summary_rows <- lapply(seq_along(models), function(i) {
  m <- models[i]
  aucs  <- get_metric(m, "auc")
  accs  <- get_metric(m, "accuracy")
  senss <- get_metric(m, "sensitivity")
  specs <- get_metric(m, "specificity")
  data.frame(
    model         = labels[i],
    mean_auc      = mean(aucs,  na.rm = TRUE),
    sd_auc        = sd(aucs,    na.rm = TRUE),
    mean_accuracy = mean(accs,  na.rm = TRUE),
    sd_accuracy   = sd(accs,    na.rm = TRUE),
    mean_sens     = mean(senss, na.rm = TRUE),
    sd_sens       = sd(senss,   na.rm = TRUE),
    mean_spec     = mean(specs, na.rm = TRUE),
    sd_spec       = sd(specs,   na.rm = TRUE)
  )
})
summary_df <- do.call(rbind, summary_rows)

message("\n=== 5-fold CV Summary ===")
print(summary_df, digits = 3)

write.csv(summary_df,    file.path(outdir, "cv_summary.csv"),      row.names = FALSE)
saveRDS(fold_results,    file.path(outdir, "cv_metrics_full.rds"))
saveRDS(list(folds = folds, y = y, samples = sids),
                         file.path(outdir, "cv_fold_assignments.rds"))

message("\nDone. Results written to: ", outdir)
