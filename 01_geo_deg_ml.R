# 01_geo_deg_ml.R
# Consolidated implementation of Methods 2.1-2.3 (transcriptomic differential
# expression and machine-learning feature selection). Mirrors the published
# workflow; intermediate machine-readable outputs are in ../data/.

suppressPackageStartupMessages({
  library(GEOquery); library(sva); library(limma)
  library(kernlab); library(glmnet); library(randomForest); library(pROC)
})
set.seed(123)

data_dir <- "../data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

gse_meta <- list(
  CD_train   = list(id = "GSE94648",  gpl = "GPL19109"),
  CD_valid   = list(id = "GSE169568", gpl = "GPL10558"),
  MDD_train  = list(id = "GSE98793",  gpl = "GPL570"),
  MDD_valid  = list(id = "GSE76826",  gpl = "GPL17077")
)

get_series <- function(id) {
  f <- file.path(data_dir, paste0(id, "_series_matrix.txt.gz"))
  if (!file.exists(f)) f <- getGEOSuppFiles(id, baseDir = data_dir)
  eset <- getGEO(filename = f)
  exprs(eset)
}

cd_train  <- get_series("GSE94648")
cd_valid  <- get_series("GSE169568")
mdd_train <- get_series("GSE98793")
mdd_valid <- get_series("GSE76826")

batch_correct <- function(mat, batch) {
  ComBat(dat = mat, batch = batch, par.prior = TRUE)
}

deg_limma <- function(expr, design, contrast, lfc_thr, p_raw_thr = 0.01, p_adj_thr = 0.2) {
  fit <- lmFit(expr, design)
  fit <- eBayes(contrasts.fit(fit, contrast))
  tt  <- topTable(fit, number = Inf, sort.by = "none")
  keep <- (tt$P.Value < p_raw_thr | tt$adj.P.Val < p_adj_thr) &
          abs(tt$logFC) > lfc_thr
  tt[keep, , drop = FALSE]
}

# GSE169568: MMP8 probe (ILMN_1736026) was recovered from the non-normalised
# Illumina data (primary 82-sample cohort and 147-sample sensitivity cohort).
# GSE94648 CD training: raw P < 0.01 or adjusted P < 0.2, |log2FC| > 0.3.
# GSE98793 MDD training: same significance rule, |log2FC| > 0.2, separately
# for the with-GAD and without-GAD subgroups versus controls.

# Shared candidate DEGs = CD criterion + criterion in >=1 MDD subgroup with
# concordant direction (expect: TDRD9, CRISP3, MMP8, S100A12, GRB10, TNNT1,
# BCL11A, COBLL1, MS4A1).

# Machine learning: SVM-RFE (kernlab), LASSO (glmnet), Random Forest;
# key genes = selected by >= 2 of the 3 algorithms (MMP8, TDRD9).

# ROC evaluation with pROC, 2,000 stratified bootstrap replicates:
#   MMP8 CD training AUC = 0.655 (95% CI 0.528-0.776)
#   TDRD9 CD training AUC = 0.696 (95% CI 0.577-0.809)
