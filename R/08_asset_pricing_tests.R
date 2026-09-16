## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026
## =============================================================================
## SECTION 8: FACTOR MODELS AND ASSET PRICING TESTS
## Baseline specification: SIC value-weighted SMQ (SMQ_sic_VW)
##
## TABLE OF CONTENTS
## -----------------
## 8.0   Setup and Data Loading
## 8.1   Model Definitions  (CAPM / FF3 / FF5 / FF5+SMQ, no orthogonalization)
## 8.2   General Helper Functions
## 8.3   Factor Summary Statistics and Spanning Tests (Tables 2 and 14)
## 8.4   Portfolio Family Specifications
## 8.5   Portfolio Helper Functions  (1D / 5x5 / 2x4x4)
## 8.6   Mean Excess Return Tables  (Tables 5 and 6)
## 8.7   Average Firm Count per Portfolio Cell
## 8.8   GRS Tests and Alpha Matrices
##         8.8.1  Per-family GRS + A|alpha| + A(R-bar^2)  (Table 4)
## 8.9   Full Loading Tables: 25-Portfolio Families  (Tables 7, and analogues)
##         One block per family x model showing alpha + all factor loadings
##         + H-L chi margin column + S-B size margin row (25-port sorts)
##         + H-L characteristic margin row (non-chi 25-port sorts)
## 8.10  Full Loading Tables: 32-Portfolio Families  (Tables 8, 9, 10, analogues)
##         One block per family x model showing alpha + all factor loadings
##         + H-L chi margin column  (chi families)
##         + H-L characteristic margin row  (chi and non-chi families)
##         Sign convention: H-L for BM, OP, chi; L-H for Inv
## 8.11  Figures
##         8.11.1 Heatmaps: FF5 alpha and loadings (Tables 14-17 style)
##         8.11.2 Line plots: alpha, HML, RMW profiles across chi
##         8.11.3 Bar charts: Small vs Big average coefficients by chi quartile
##         8.11.4 Bar chart: H-L chi spreads by characteristic quartile
##         8.11.5 Bar chart: Average H-L chi spread (compact)
## =============================================================================

### 8.0 Setup and Data Loading ################################################
## Clear Everything
cat("\014")
rm(list = ls())
graphics.off()

## Project paths
source(file.path("R", "00_config.R"))

## Download Libaries
library(data.table)
library(lubridate)
library(lmtest)
library(sandwich)
library(ggplot2)

## Sample Period and Lags
START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")
NW_LAGS <- 12 # Newey-West lags throughout

## Load master panel (factors + test assets merged in Section 5)
panel <- readRDS("master_factors_testassets.rds")
setDT(panel)
panel[, mdate := as.Date(mdate)]

factor_test_panel <- panel[mdate >= START_DATE & mdate <= END_DATE]
rf_monthly <- unique(factor_test_panel[, .(mdate, RF)])

## Load labelled portfolio data (for firm counts)
labelled_data <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(labelled_data)
labelled_data[, mdate := as.Date(mdate)]
labelled_data <- labelled_data[mdate >= START_DATE & mdate <= END_DATE]

## Output directory for figures
fig_dir <- file.path(FIGURES_DIR, "section8")
if (!dir.exists(fig_dir)) dir.create(fig_dir)

save_fig <- function(p, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
  invisible(p)
}


### 8.1 Model Definitions ######################################################
## Four models, NO orthogonalization of any factor.
## FF5 uses raw HML. FF5+SMQ adds SMQ_sic_VW directly on top of FF5.

models <- list(
  CAPM    = c("MKT"),
  FF3     = c("MKT", "SMB", "HML"),
  FF5     = c("MKT", "SMB", "HML", "RMW", "CMA"),
  FF5_SMQ = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_VW")
)

### 8.2 General Helper Functions ##############################################
## Star labels from p-values
star_p <- function(p) {
  fifelse(
    p < 0.01, "***",
    fifelse(
      p < 0.05, "**",
      fifelse(p < 0.10, "*", "")
    )
  )
}

## Star labels from t-statistics
star_t <- function(t) {
  if (is.na(t)) {
    return("")
  }
  if (abs(t) >= 2.58) {
    return("***")
  }
  if (abs(t) >= 1.96) {
    return("**")
  }
  if (abs(t) >= 1.65) {
    return("*")
  }
  return("")
}

## Round all numeric columns in a data.table
round_numeric_dt <- function(dt, digits = 4) {
  out <- copy(dt)
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  out
}

## Detect portfolio columns by regex, with optional count check
get_cols <- function(x, pattern, expected_n = NULL, label = pattern) {
  cols <- grep(pattern, x, value = TRUE)
  if (length(cols) == 0) stop("No columns matched: ", label)
  if (!is.null(expected_n) && length(cols) != expected_n) {
    stop("Expected ", expected_n, " cols for ", label, ", got ", length(cols))
  }
  cols
}

## Time-series OLS with Newey-West SEs (prewhite = FALSE, adjust = TRUE)
## Returns a data.table with term, estimate, tstat, pval, adj_r2, stars
run_ts_reg_nw <- function(y, X, lags = NW_LAGS) {
  dt <- data.table(y = y)
  dt <- cbind(dt, as.data.table(X))
  dt <- dt[complete.cases(dt)]
  if (nrow(dt) < 60) {
    return(NULL)
  }
  fml <- as.formula(paste("y ~", paste(names(X), collapse = " + ")))
  fit <- lm(fml, data = dt)
  nw <- NeweyWest(fit, lag = lags, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov. = nw)
  out <- data.table(
    term     = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    tstat    = as.numeric(ct[, 3]),
    pval     = as.numeric(ct[, 4]),
    adj_r2   = summary(fit)$adj.r.squared
  )
  out[, stars := star_p(pval)]
  out[]
}

## Newey-West t-stat for a mean (regress on intercept only)
nw_tstat_mean <- function(x, lags = NW_LAGS) {
  x <- x[!is.na(x)]
  if (length(x) < 60) {
    return(NA_real_)
  }
  fit <- lm(x ~ 1)
  nw <- NeweyWest(fit, lag = lags, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov. = nw)
  as.numeric(ct[1, "t value"])
}

## Format a coefficient cell: "est***\n(tstat)"
fmt_cell <- function(est, tstat, stars,
                     scale_est = 1, d_est = 2, d_t = 2) {
  out <- sprintf(
    paste0("%.", d_est, "f%s\n(%.", d_t, "f)"),
    scale_est * est,
    stars,
    tstat
  )
  out[is.na(est) | is.na(tstat)] <- ""
  out
}

## GRS F-test for a block of test assets
## Returns data.table with GRS, GRS_pval, mean_abs_alpha, avg_adj_r2
run_grs <- function(dt, asset_cols, factor_cols) {
  tmp <- copy(dt[, c("RF", asset_cols, factor_cols), with = FALSE])
  for (col in asset_cols) tmp[, (col) := get(col) - RF]
  tmp[, RF := NULL]
  tmp <- tmp[complete.cases(tmp)]
  Tn <- nrow(tmp)
  N <- length(asset_cols)
  K <- length(factor_cols)
  if (Tn <= N + K) {
    return(NULL)
  }
  R <- as.matrix(tmp[, ..asset_cols])
  F <- as.matrix(tmp[, ..factor_cols])
  X <- cbind(1, F)
  B <- solve(crossprod(X)) %*% crossprod(X, R)
  alpha <- as.numeric(B[1, ])
  E <- R - X %*% B
  Sigma_e <- crossprod(E) / (Tn - K - 1)
  Sigma_f <- cov(F)
  mu_f <- colMeans(F)
  if (qr(Sigma_e)$rank < N || qr(Sigma_f)$rank < K) {
    return(NULL)
  }
  grs_num <- ((Tn - N - K) / N) *
    t(matrix(alpha)) %*% solve(Sigma_e) %*% matrix(alpha)
  grs_den <- 1 + t(matrix(mu_f)) %*% solve(Sigma_f) %*% matrix(mu_f)
  GRS <- as.numeric(grs_num / grs_den)
  GRS_pval <- 1 - pf(GRS, N, Tn - N - K)
  # Average adjusted R2
  fitted_vals <- X %*% B
  resid_vals <- R - fitted_vals
  y_bar <- matrix(colMeans(R), nrow = Tn, ncol = N, byrow = TRUE)
  tss <- colSums((R - y_bar)^2)
  rss <- colSums(resid_vals^2)
  adj_r2 <- 1 - (rss / (Tn - K - 1)) / (tss / (Tn - 1))
  mean_excess <- colMeans(R)
  data.table(
    T = Tn, N = N, K = K,
    GRS = GRS,
    GRS_pval = GRS_pval,
    mean_abs_alpha = 100 * mean(abs(alpha)),
    avg_adj_r2 = mean(adj_r2, na.rm = TRUE)
  )
}

## Run NW regression for a single long-short series column
run_ls_reg_nw <- function(panel, col, rhs, lags = NW_LAGS) {
  y <- panel[[col]]
  if (is.null(y)) {
    return(NULL)
  }
  X <- panel[, ..rhs]
  run_ts_reg_nw(y, X, lags)
}

## Format a single alpha cell from a regression result
fmt_alpha_cell <- function(reg, d_est = 2, d_t = 2) {
  if (is.null(reg)) {
    return("")
  }
  a <- reg[term == "(Intercept)"]
  if (nrow(a) == 0) {
    return("")
  }
  fmt_cell(a$estimate, a$tstat, a$stars,
    scale_est = 100,
    d_est = d_est, d_t = d_t
  )
}


### 8.3 Factor Summary Statistics and Spanning Tests ##########################
## Produces:
##   Table 2  -- Panel A: mean, SD, Sharpe for all 6 baseline factors
##              Panel B: factor-on-factor spanning regressions (all 6 factors)
##   Table 14 -- FF5-only spanning regressions (5 factors)

## Panel A: summary statistics for MKT, SMB, HML, RMW, CMA, SMQ_sic_VW
factor_set_6 <- c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_VW")
factor_set_5 <- c("MKT", "SMB", "HML", "RMW", "CMA")

make_factor_summary <- function(panel, factor_vec) {
  tab <- data.table(Statistic = c("Mean excess return", "Std. Dev.", "Sharpe ratio"))
  for (f in factor_vec) {
    x <- panel[[f]]
    mu <- 100 * mean(x, na.rm = TRUE)
    sd <- 100 * sd(x, na.rm = TRUE)
    sr <- sqrt(12) * mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)
    ts <- nw_tstat_mean(x, NW_LAGS)
    st <- star_t(ts)
    tab[, (f) := c(
      sprintf("%.2f%s\n(%.2f)", mu, st, ts),
      sprintf("%.2f", sd),
      sprintf("%.2f", sr)
    )]
  }
  tab
}

## Helper: run factor-on-factor spanning regressions
## Each factor is regressed on all remaining factors in the set
make_spanning_table <- function(panel, factor_vec) {
  res <- list()
  for (y_name in factor_vec) {
    x_names <- setdiff(factor_vec, y_name)
    tmp <- panel[, c(y_name, x_names), with = FALSE]
    tmp <- tmp[complete.cases(tmp)]
    fml <- as.formula(paste(y_name, "~", paste(x_names, collapse = " + ")))
    fit <- lm(fml, data = tmp)
    nw <- NeweyWest(fit, lag = NW_LAGS, prewhite = FALSE, adjust = TRUE)
    ct <- coeftest(fit, vcov. = nw)
    row <- data.table(dep = y_name)
    # Intercept (alpha)
    a <- ct["(Intercept)", ]
    row[, alpha := fmt_cell(a[1], a[3], star_p(a[4]), scale_est = 100)]
    # Slope coefficients (reported as standard-beta-unit loadings)
    for (xx in x_names) {
      b <- ct[xx, ]
      row[, (xx) := fmt_cell(b[1], b[3], star_p(b[4]))]
    }
    row[, Adj_R2 := sprintf("%.2f", summary(fit)$adj.r.squared)]
    res[[y_name]] <- row
  }
  rbindlist(res, fill = TRUE)
}


## Correlations between the six baseline factors.
## This is descriptive only and does not replace the spanning regressions.
make_factor_corr_table <- function(panel, factor_vec) {
  tmp <- panel[, ..factor_vec]
  tmp <- tmp[complete.cases(tmp)]

  corr_mat <- cor(tmp)

  ## Round for printed output
  corr_tab <- as.data.table(round(corr_mat, 2), keep.rownames = "Factor")

  list(
    corr_mat = corr_mat,
    corr_tab = corr_tab
  )
}

## Summary Statistics
table2_panel_A <- make_factor_summary(factor_test_panel, factor_set_6)

## Factor Spanning Tests for FF5 + SMQ
table2_panel_B <- make_spanning_table(factor_test_panel, factor_set_6)

## Summary Statistics and Spanning Tests
table14_ff5_spanning <- make_spanning_table(factor_test_panel, factor_set_5)

## Correlation Matrix
factor_corr <- make_factor_corr_table(factor_test_panel, factor_set_6)

cat("\n============================================================\n")
cat("FACTOR SUMMARY STATISTICS\n")
cat("============================================================\n")
print(table2_panel_A)

cat("\n============================================================\n")
cat("FACTOR SPANNING REGRESSIONS (FF5 + SMQ FACTOR SET)\n")
cat("============================================================\n")
print(table2_panel_B)

cat("\n============================================================\n")
cat("FACTOR SPANNING REGRESSIONS (FF5 FACTOR SET)\n")
cat("============================================================\n")
print(table14_ff5_spanning)

cat("\n============================================================\n")
cat("FACTOR CORRELATION MATRIX\n")
cat("============================================================\n")
print(factor_corr$corr_tab)

saveRDS(table2_panel_A, "out_table2_panel_A.rds")
saveRDS(table2_panel_B, "out_table2_panel_B.rds")
saveRDS(table14_ff5_spanning, "out_table14_ff5_spanning.rds")
saveRDS(factor_corr$corr_mat, "out_factor_correlation_matrix.rds")
saveRDS(factor_corr$corr_tab, "out_factor_correlation_table.rds")

### 8.4 Portfolio Family Specifications #######################################
## We work with two groups of test assets:
##   (a) 25-portfolio families (5x5 sorts): BM, OP, INV, CHI_SIC,
##                                          BM_CHI_SIC, OP_CHI_SIC, INV_CHI_SIC
##   (b) 32-portfolio families (2x4x4 sorts): BM_OP, BM_INV, OP_INV,
##                                             OP_CHI_SIC, BM_CHI_SIC, INV_CHI_SIC
##
## Each 5x5 entry specifies the column prefix and row/col codes used in
## factor_test_panel column names (e.g. "BM_S3V4").
## Each 2x4x4 entry specifies the column prefix (e.g. "BM_OP_S32").

families_25 <- list(
  BM          = list(prefix = "BM_", row = "S", col = "V"),
  OP          = list(prefix = "OP_", row = "S", col = "P"),
  INV         = list(prefix = "INV_", row = "S", col = "I"),
  CHI_SIC     = list(prefix = "CHI_SIC_", row = "S", col = "C"),
  BM_CHI_SIC  = list(prefix = "BM_CHI_SIC_", row = "B", col = "C"),
  OP_CHI_SIC  = list(prefix = "OP_CHI_SIC_", row = "O", col = "C"),
  INV_CHI_SIC = list(prefix = "INV_CHI_SIC_", row = "I", col = "C")
)

families_32 <- list(
  BM_OP       = list(prefix = "BM_OP_"),
  BM_INV      = list(prefix = "BM_INV_"),
  OP_INV      = list(prefix = "OP_INV_"),
  OP_CHI_SIC  = list(prefix = "OP_CHI_SIC_"),
  BM_CHI_SIC  = list(prefix = "BM_CHI_SIC_2X4X4_"),
  INV_CHI_SIC = list(prefix = "INV_CHI_SIC_2X4X4_")
)

## Families that include chi as the second sort dimension (column = chi)
chi_families_25 <- c("CHI_SIC", "BM_CHI_SIC", "OP_CHI_SIC", "INV_CHI_SIC")
chi_families_32 <- c("OP_CHI_SIC", "BM_CHI_SIC", "INV_CHI_SIC")

## Inv families where the characteristic spread is L-H (not H-L)
inv_families <- c("INV", "INV_CHI_SIC", "BM_INV", "OP_INV")


### 8.5 Portfolio Helper Functions ############################################

## ---- 5x5 helpers -----------------------------------------------------------

## Convert 5x5 wide panel columns to long format
wide_5x5_to_long <- function(panel, prefix, row_code, col_code) {
  cols <- get_cols(
    names(panel),
    paste0("^", prefix, row_code, "[1-5]", col_code, "[1-5]$"),
    25, paste0(prefix, " 5x5")
  )
  long <- melt(panel[, c("mdate", cols), with = FALSE],
    id.vars = "mdate", variable.name = "portfolio", value.name = "vwret"
  )
  long <- long[!is.na(vwret)]
  x <- sub(paste0("^", prefix), "", as.character(long$portfolio))
  long[, row5 := substr(x, 1, 2)] # e.g. "S3"
  long[, col5 := substr(x, 3, 4)] # e.g. "V4"
  long[]
}

## Mean excess return 5x5 matrix
make_5x5_mean_excess <- function(long_dt, rf_dt, row_levels, col_levels) {
  dt <- merge(long_dt, rf_dt, by = "mdate", all.x = TRUE)
  dt[, excess := vwret - RF]
  m <- dt[, .(val = 100 * mean(excess, na.rm = TRUE)), by = .(row5, col5)]
  m[, row5 := factor(row5, levels = row_levels)]
  m[, col5 := factor(col5, levels = col_levels)]
  mat <- dcast(m, row5 ~ col5, value.var = "val")
  out <- as.matrix(mat[, -1, with = FALSE])
  rownames(out) <- mat$row5
  out
}

## Run OLS-NW regressions for all 25 portfolios in one family
run_regs_5x5 <- function(panel, prefix, row_code, col_code, rhs) {
  cols <- get_cols(
    names(panel),
    paste0("^", prefix, row_code, "[1-5]", col_code, "[1-5]$"),
    25, prefix
  )
  X <- panel[, ..rhs]
  res <- lapply(cols, function(col) {
    pat <- sub(paste0("^", prefix), "", col)
    r5 <- substr(pat, 1, 2)
    c5 <- substr(pat, 3, 4)
    y <- panel[[col]] - panel[["RF"]]
    est <- run_ts_reg_nw(y, X)
    if (is.null(est)) {
      return(NULL)
    }
    est[, `:=`(portfolio = col, row5 = r5, col5 = c5)]
    est
  })
  rbindlist(res, fill = TRUE)
}

## Extract a 5x5 coefficient matrix from regression output
term_to_5x5 <- function(reg_dt, term_name, row_levels, col_levels,
                        scale_est = 1, d_est = 2, d_t = 2) {
  tmp <- reg_dt[term == term_name, .(row5, col5, estimate, tstat, stars)]
  if (nrow(tmp) == 0) {
    return(matrix("", 5, 5))
  }
  tmp[, cell := fmt_cell(estimate, tstat, stars, scale_est, d_est, d_t)]
  tmp[, row5 := factor(row5, levels = row_levels)]
  tmp[, col5 := factor(col5, levels = col_levels)]
  mat <- dcast(tmp, row5 ~ col5, value.var = "cell", fill = "")
  out <- as.matrix(mat[, -1, with = FALSE])
  rownames(out) <- mat$row5
  out
}

## ---- 2x4x4 helpers ---------------------------------------------------------

## Convert 2x4x4 wide panel columns to long format
wide_2x4x4_to_long <- function(panel, prefix) {
  cols <- get_cols(
    names(panel),
    paste0("^", prefix, "[SB][1-4][1-4]$"),
    32, paste0(prefix, " 2x4x4")
  )
  long <- melt(panel[, c("mdate", cols), with = FALSE],
    id.vars = "mdate", variable.name = "portfolio", value.name = "vwret"
  )
  long <- long[!is.na(vwret)]
  x <- sub(paste0("^", prefix), "", as.character(long$portfolio))
  long[, size := substr(x, 1, 1)] # "S" or "B"
  long[, c1 := substr(x, 2, 2)] # row characteristic quartile
  long[, c2 := substr(x, 3, 3)] # column characteristic quartile
  long[]
}

## Mean excess return 2x4x4 matrix (returns list with S and B)
make_2x4x4_mean_excess <- function(long_dt, rf_dt) {
  dt <- merge(long_dt, rf_dt, by = "mdate", all.x = TRUE)
  dt[, excess := vwret - RF]
  m <- dt[, .(val = 100 * mean(excess, na.rm = TRUE)), by = .(size, c1, c2)]
  q_levels <- as.character(1:4)
  out <- list()
  for (sz in c("S", "B")) {
    sub <- m[size == sz]
    sub[, c1 := factor(c1, levels = q_levels)]
    sub[, c2 := factor(c2, levels = q_levels)]
    mat <- dcast(sub, c1 ~ c2, value.var = "val", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- paste0("Q", q_levels)
    colnames(mat_out) <- paste0("Q", q_levels)
    out[[sz]] <- mat_out
  }
  out
}

## Run OLS-NW regressions for all 32 portfolios in one family
run_regs_2x4x4 <- function(panel, prefix, rhs) {
  cols <- get_cols(
    names(panel),
    paste0("^", prefix, "[SB][1-4][1-4]$"),
    32, prefix
  )
  X <- panel[, ..rhs]
  res <- lapply(cols, function(col) {
    pat <- sub(paste0("^", prefix), "", col)
    size <- substr(pat, 1, 1)
    c1 <- substr(pat, 2, 2)
    c2 <- substr(pat, 3, 3)
    y <- panel[[col]] - panel[["RF"]]
    est <- run_ts_reg_nw(y, X)
    if (is.null(est)) {
      return(NULL)
    }
    est[, `:=`(portfolio = col, size = size, c1 = c1, c2 = c2)]
    est
  })
  rbindlist(res, fill = TRUE)
}

## Extract a 4x4 coefficient matrix for one size group
term_to_4x4 <- function(reg_dt, term_name, size_flag,
                        scale_est = 1, d_est = 2, d_t = 2) {
  tmp <- reg_dt[
    term == term_name & size == size_flag,
    .(c1, c2, estimate, tstat, stars)
  ]
  if (nrow(tmp) == 0) {
    return(matrix("", 4, 4))
  }
  tmp[, cell := fmt_cell(estimate, tstat, stars, scale_est, d_est, d_t)]
  tmp[, c1 := factor(c1, levels = as.character(1:4))]
  tmp[, c2 := factor(c2, levels = as.character(1:4))]
  mat <- dcast(tmp, c1 ~ c2, value.var = "cell", fill = "")
  out <- as.matrix(mat[, -1, with = FALSE])
  rownames(out) <- paste0("Q", mat$c1)
  out
}


### 8.6 Mean Excess Return Tables (Tables 5 and 6) ###########################
## Table 5: mean monthly excess returns for each 25-portfolio family
## Table 6: mean monthly excess returns for each 32-portfolio family

cat("\n============================================================\n")
cat("TABLE 5: MEAN MONTHLY EXCESS RETURNS - 25-PORTFOLIO FAMILIES\n")
cat("============================================================\n")

mean_excess_25 <- list()

for (fam in names(families_25)) {
  spec <- families_25[[fam]]
  long <- wide_5x5_to_long(factor_test_panel, spec$prefix, spec$row, spec$col)

  row_levels <- paste0(spec$row, 1:5)
  col_levels <- paste0(spec$col, 1:5)

  mat <- make_5x5_mean_excess(long, rf_monthly, row_levels, col_levels)
  mean_excess_25[[fam]] <- mat

  cat("\n---", fam, "---\n")
  print(round(mat, 2))
}

cat("\n============================================================\n")
cat("TABLE 6: MEAN MONTHLY EXCESS RETURNS - 32-PORTFOLIO FAMILIES\n")
cat("============================================================\n")

mean_excess_32 <- list()

for (fam in names(families_32)) {
  spec <- families_32[[fam]]
  long <- wide_2x4x4_to_long(factor_test_panel, spec$prefix)

  res <- make_2x4x4_mean_excess(long, rf_monthly)
  mean_excess_32[[fam]] <- res

  cat("\n---", fam, "---\n")
  cat("Small\n")
  print(round(res$S, 2))
  cat("Big\n")
  print(round(res$B, 2))
}

saveRDS(mean_excess_25, "out_table5_mean_excess_25port.rds")
saveRDS(mean_excess_32, "out_table6_mean_excess_32port.rds")


### 8.7 Average Number of Firms per Portfolio Cell ############################
## Helps assess portfolio diversification.
## For each family, we compute the time-averaged number of unique firms
## in each cell over the sample period.

cat("\n============================================================\n")
cat("AVERAGE NUMBER OF FIRMS PER PORTFOLIO CELL\n")
cat("============================================================\n")

## ---- 25-portfolio families -------------------------------------------------
## Portfolio labels in labelled_data use the same naming conventions.
## We use the chi1d (5-bin chi) and the standard size/bm/op/inv 5-bin sorts.

firm_count_25 <- list()

## Map from family name to the labelled_data column and its unique values
family_label_cols_25 <- list(
  BM          = list(row_col = c("size5", "value5"), row_pfx = "S", col_pfx = "V"),
  OP          = list(row_col = c("size5", "prof5"), row_pfx = "S", col_pfx = "P"),
  INV         = list(row_col = c("size5", "invest5"), row_pfx = "S", col_pfx = "I"),
  CHI_SIC     = list(row_col = c("size5", "chi5"), row_pfx = "S", col_pfx = "C"),
  BM_CHI_SIC  = list(row_col = c("bm5_chi", "chi5"), row_pfx = "B", col_pfx = "C"),
  OP_CHI_SIC  = list(row_col = c("op5_chi", "chi5"), row_pfx = "O", col_pfx = "C"),
  INV_CHI_SIC = list(row_col = c("inv5_chi", "chi5"), row_pfx = "I", col_pfx = "C")
)

for (fam in names(family_label_cols_25)) {
  cfg <- family_label_cols_25[[fam]]
  rc <- cfg$row_col

  if (!all(rc %in% names(labelled_data))) {
    cat("  Skipping", fam, "(label columns not found in labelled_data)\n")
    next
  }

  # Average unique firms per cell per month, then average across months
  counts <- labelled_data[
    !is.na(get(rc[1])) & !is.na(get(rc[2])),
    .(n_firms = uniqueN(permno)),
    by = c("mdate", rc)
  ][,
    .(avg_n = mean(n_firms, na.rm = TRUE)),
    by = rc
  ]

  setnames(counts, rc, c("row_dim", "col_dim"))

  row_levels <- paste0(cfg$row_pfx, 1:5)
  col_levels <- paste0(cfg$col_pfx, 1:5)

  counts[, row_dim := factor(row_dim, levels = row_levels)]
  counts[, col_dim := factor(col_dim, levels = col_levels)]

  mat <- dcast(counts, row_dim ~ col_dim, value.var = "avg_n", drop = FALSE)
  mat_out <- round(as.matrix(mat[, -1, with = FALSE]), 0)
  rownames(mat_out) <- as.character(mat$row_dim)

  firm_count_25[[fam]] <- mat_out

  cat("\n---", fam, "---\n")
  print(mat_out)
}

## ---- 32-portfolio families -------------------------------------------------
## 2x4x4 portfolio labels are stored as port_<TAG>_2x4x4 columns
firm_count_32 <- list()

port_col_map_32 <- list(
  BM_OP       = "port_BM_OP_2x4x4",
  BM_INV      = "port_BM_INV_2x4x4",
  OP_INV      = "port_OP_INV_2x4x4",
  OP_CHI_SIC  = "port_OP_CHI_SIC_2x4x4",
  BM_CHI_SIC  = "port_BM_CHI_SIC_2x4x4",
  INV_CHI_SIC = "port_INV_CHI_SIC_2x4x4"
)

for (fam in names(port_col_map_32)) {
  pc <- port_col_map_32[[fam]]

  if (!pc %in% names(labelled_data)) {
    cat("  Skipping", fam, "(port column not found)\n")
    next
  }

  counts <- labelled_data[
    !is.na(get(pc)),
    .(n_firms = uniqueN(permno)),
    by = c("mdate", pc)
  ][,
    .(avg_n = mean(n_firms, na.rm = TRUE)),
    by = pc
  ]

  counts[, size := substr(get(pc), 1, 1)]
  counts[, c1 := substr(get(pc), 2, 2)]
  counts[, c2 := substr(get(pc), 3, 3)]

  q_levels <- as.character(1:4)
  out <- list()
  for (sz in c("S", "B")) {
    sub <- counts[size == sz]
    sub[, c1 := factor(c1, levels = q_levels)]
    sub[, c2 := factor(c2, levels = q_levels)]
    mat <- dcast(sub, c1 ~ c2, value.var = "avg_n", drop = FALSE)
    mat_out <- round(as.matrix(mat[, -1, with = FALSE]), 0)
    rownames(mat_out) <- paste0("Q", q_levels)
    colnames(mat_out) <- paste0("Q", q_levels)
    out[[sz]] <- mat_out
  }

  firm_count_32[[fam]] <- out

  cat("\n---", fam, "---\n")
  cat("Small\n")
  print(out$S)
  cat("Big\n")
  print(out$B)
}

saveRDS(firm_count_25, "out_avg_firm_count_25port.rds")
saveRDS(firm_count_32, "out_avg_firm_count_32port.rds")


### 8.8 GRS Tests (Table 4) ###################################################
## For each portfolio family, and for each model, we report:
##   GRS statistic and p-value
##   Average absolute alpha A|alpha_i|
##   Average adjusted R^2  A(R-bar^2)

## Collect all 25-portfolio asset column sets
all_25_cols <- c(
  get_cols(names(factor_test_panel), "^BM_S[1-5]V[1-5]$", 25, "BM"),
  get_cols(names(factor_test_panel), "^OP_S[1-5]P[1-5]$", 25, "OP"),
  get_cols(names(factor_test_panel), "^INV_S[1-5]I[1-5]$", 25, "INV"),
  get_cols(names(factor_test_panel), "^CHI_SIC_S[1-5]C[1-5]$", 25, "CHI_SIC"),
  get_cols(names(factor_test_panel), "^BM_CHI_SIC_B[1-5]C[1-5]$", 25, "BM_CHI_SIC"),
  get_cols(names(factor_test_panel), "^OP_CHI_SIC_O[1-5]C[1-5]$", 25, "OP_CHI_SIC"),
  get_cols(names(factor_test_panel), "^INV_CHI_SIC_I[1-5]C[1-5]$", 25, "INV_CHI_SIC")
)

all_32_cols <- c(
  get_cols(names(factor_test_panel), "^BM_OP_[SB][1-4][1-4]$", 32, "BM_OP"),
  get_cols(names(factor_test_panel), "^BM_INV_[SB][1-4][1-4]$", 32, "BM_INV"),
  get_cols(names(factor_test_panel), "^OP_INV_[SB][1-4][1-4]$", 32, "OP_INV"),
  get_cols(names(factor_test_panel), "^OP_CHI_SIC_[SB][1-4][1-4]$", 32, "OP_CHI_SIC_32"),
  get_cols(names(factor_test_panel), "^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", 32, "BM_CHI_SIC_32"),
  get_cols(names(factor_test_panel), "^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", 32, "INV_CHI_SIC_32")
)

## ---- 8.8.1 Per-family GRS (Table 4) ----------------------------------------
grs_results <- list()

## 25-portfolio families
asset_col_map_25 <- list(
  BM          = get_cols(names(factor_test_panel), "^BM_S[1-5]V[1-5]$", 25),
  OP          = get_cols(names(factor_test_panel), "^OP_S[1-5]P[1-5]$", 25),
  INV         = get_cols(names(factor_test_panel), "^INV_S[1-5]I[1-5]$", 25),
  CHI_SIC     = get_cols(names(factor_test_panel), "^CHI_SIC_S[1-5]C[1-5]$", 25),
  BM_CHI_SIC  = get_cols(names(factor_test_panel), "^BM_CHI_SIC_B[1-5]C[1-5]$", 25),
  OP_CHI_SIC  = get_cols(names(factor_test_panel), "^OP_CHI_SIC_O[1-5]C[1-5]$", 25),
  INV_CHI_SIC = get_cols(names(factor_test_panel), "^INV_CHI_SIC_I[1-5]C[1-5]$", 25)
)

asset_col_map_32 <- list(
  BM_OP       = get_cols(names(factor_test_panel), "^BM_OP_[SB][1-4][1-4]$", 32),
  BM_INV      = get_cols(names(factor_test_panel), "^BM_INV_[SB][1-4][1-4]$", 32),
  OP_INV      = get_cols(names(factor_test_panel), "^OP_INV_[SB][1-4][1-4]$", 32),
  OP_CHI_SIC  = get_cols(names(factor_test_panel), "^OP_CHI_SIC_[SB][1-4][1-4]$", 32),
  BM_CHI_SIC  = get_cols(names(factor_test_panel), "^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", 32),
  INV_CHI_SIC = get_cols(names(factor_test_panel), "^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", 32)
)

## Run GRS for each family x model combination
run_grs_family <- function(asset_col_map, family_label, n_port) {
  results <- list()
  for (fam in names(asset_col_map)) {
    for (mod in names(models)) {
      res <- run_grs(factor_test_panel, asset_col_map[[fam]], models[[mod]])
      if (!is.null(res)) {
        res[, `:=`(
          family = fam, model = mod,
          port_label = family_label, n_portfolios = n_port
        )]
        results[[paste(fam, mod)]] <- res
      }
    }
  }
  rbindlist(results, fill = TRUE)
}

grs_25 <- run_grs_family(asset_col_map_25, "25-portfolio", 25)
grs_32 <- run_grs_family(asset_col_map_32, "32-portfolio", 32)

## Display for 25 portfolios
cat("\n============================================================\n")
cat("TABLE 4: GRS TESTS - 25-PORTFOLIO FAMILIES\n")
cat("Columns: GRS | p-val | A|alpha| (%) | A(R-bar^2)\n")
cat("============================================================\n")

model_order <- c("CAPM", "FF3", "FF5", "FF5_SMQ")

for (fam in names(asset_col_map_25)) {
  cat("\n---", fam, "---\n")

  sub <- copy(grs_25[family == fam & model %in% model_order])
  sub[, model := factor(model, levels = model_order)]
  setorder(sub, model)

  sub[, GRS_star := paste0(sprintf("%.2f", GRS), star_p(GRS_pval))]
  sub[, GRS_pval_fmt := sprintf("%.3f", GRS_pval)]
  sub[, mean_abs_alpha_fmt := sprintf("%.2f", mean_abs_alpha)]
  sub[, avg_adj_r2_fmt := sprintf("%.2f", avg_adj_r2)]

  print(sub[, .(
    model,
    GRS = GRS_star,
    GRS_pval = GRS_pval_fmt,
    mean_abs_alpha = mean_abs_alpha_fmt,
    avg_adj_r2 = avg_adj_r2_fmt
  )])
}

## Display for 32 portfolios
cat("\n============================================================\n")
cat("TABLE 4: GRS TESTS - 32-PORTFOLIO FAMILIES\n")
cat("Columns: GRS | p-val | A|alpha| (%) | A(R-bar^2)\n")
cat("============================================================\n")

for (fam in names(asset_col_map_32)) {
  cat("\n---", fam, "---\n")

  sub <- copy(grs_32[family == fam & model %in% model_order])
  sub[, model := factor(model, levels = model_order)]
  setorder(sub, model)

  sub[, GRS_star := paste0(sprintf("%.2f", GRS), star_p(GRS_pval))]
  sub[, GRS_pval_fmt := sprintf("%.3f", GRS_pval)]
  sub[, mean_abs_alpha_fmt := sprintf("%.2f", mean_abs_alpha)]
  sub[, avg_adj_r2_fmt := sprintf("%.2f", avg_adj_r2)]

  print(sub[, .(
    model,
    GRS = GRS_star,
    GRS_pval = GRS_pval_fmt,
    mean_abs_alpha = mean_abs_alpha_fmt,
    avg_adj_r2 = avg_adj_r2_fmt
  )])
}

## Save outputs
saveRDS(grs_25, "out_table4_grs_25port.rds")
saveRDS(grs_32, "out_table4_grs_32port.rds")


### 8.9 25-Portfolio Families: Alphas and loadings ############################
## Structure mirrors the paper's Tables 7-10 for 25-port sorts.
## For each family x model we display panels for:
##   - Alpha (scaled to %, NW t-stats below)
##   - Each factor loading (MKT for CAPM; MKT+SMB+HML for FF3; etc.)
##   - Adj R^2
##
## Margin additions:
##   - Rightmost column:
##       H-L chi  for chi-column families
##       H-L char for BM and OP column families
##       L-H char for Inv column families
##
##   - Bottom row:
##       S-B      for size-row families
##       H-L char for BM and OP row families
##       L-H char for Inv row families
##
## Sign conventions:
##   - H-L for BM, OP, chi
##   - L-H for Inv
##   - S-B for size
##
## The bottom-right corner is left blank because this would be a
## double-difference portfolio, which is not constructed here.


### 8.9.1 Long-short helpers ##################################################

## Build one long-short return series directly from two portfolio columns.
## The input columns should already be portfolio returns, not excess returns.
## Since both sides are portfolio returns, RF cancels in the long-short spread.
compute_ls_series_5x5 <- function(panel, long_col, short_col) {
  if (!long_col %in% names(panel) || !short_col %in% names(panel)) {
    return(NULL)
  }

  data.table(
    mdate = panel$mdate,
    ret   = panel[[long_col]] - panel[[short_col]]
  )
}

## Run NW regression on an already constructed long-short return series
run_ls_series_nw <- function(series_dt, panel, rhs, lags = NW_LAGS) {
  if (is.null(series_dt)) {
    return(NULL)
  }

  dt <- merge(
    series_dt,
    panel[, c("mdate", rhs), with = FALSE],
    by = "mdate",
    all.x = TRUE
  )

  y <- dt$ret
  X <- dt[, ..rhs]

  run_ts_reg_nw(y, X, lags)
}

## Format any term from a long-short regression
fmt_term_cell <- function(reg, term_name, scale_est = 1, d_est = 2, d_t = 2) {
  if (is.null(reg)) {
    return("")
  }

  b <- reg[term == term_name]
  if (nrow(b) == 0) {
    return("")
  }

  fmt_cell(
    est       = b$estimate,
    tstat     = b$tstat,
    stars     = b$stars,
    scale_est = scale_est,
    d_est     = d_est,
    d_t       = d_t
  )
}


### 8.9.2 Helper: build full loading display for one 5x5 family x model #######

make_full_loading_table_25 <- function(panel, fam, spec, model_name, rhs) {
  row_levels <- paste0(spec$row, 1:5)
  col_levels <- paste0(spec$col, 1:5)

  reg <- run_regs_5x5(panel, spec$prefix, spec$row, spec$col, rhs)
  if (is.null(reg) || nrow(reg) == 0) {
    return(NULL)
  }

  ## Terms to display
  all_terms <- c("(Intercept)", rhs)

  ## Alpha is reported in monthly %, factor loadings are not scaled
  scale_for <- function(tt) ifelse(tt == "(Intercept)", 100, 1)

  ## Rightmost column label
  ## This is the spread across the column dimension within each row.
  col_spread_label <- if (spec$col == "C") {
    "H-L chi"
  } else if (spec$col == "I") {
    "L-H char"
  } else {
    "H-L char"
  }

  ## Bottom row label
  ## This is the spread across the row dimension within each column.
  row_spread_label <- if (spec$row == "S") {
    "S-B"
  } else if (spec$row == "I") {
    "L-H char"
  } else {
    "H-L char"
  }

  result <- list()

  for (tt in all_terms) {
    ## Base 5x5 coefficient matrix
    base_mat <- term_to_5x5(
      reg,
      term_name  = tt,
      row_levels = row_levels,
      col_levels = col_levels,
      scale_est  = scale_for(tt)
    )

    ## ------------------------------------------------------------
    ## Rightmost column: spread across columns within each row
    ##
    ## Examples:
    ##   BM:          V5 - V1 within each size row
    ##   OP:          P5 - P1 within each size row
    ##   INV:         I1 - I5 within each size row
    ##   CHI_SIC:     C5 - C1 within each size row
    ##   BM_CHI_SIC:  C5 - C1 within each B/M row
    ##   OP_CHI_SIC:  C5 - C1 within each OP row
    ##   INV_CHI_SIC: C5 - C1 within each Inv row
    ## ------------------------------------------------------------

    col_spread <- sapply(row_levels, function(rr) {
      low_col <- paste0(spec$prefix, rr, col_levels[1])
      high_col <- paste0(spec$prefix, rr, col_levels[5])

      if (spec$col == "I") {
        ## Investment effect is low-minus-high
        series_dt <- compute_ls_series_5x5(panel, low_col, high_col)
      } else {
        ## BM, OP and chi use high-minus-low
        series_dt <- compute_ls_series_5x5(panel, high_col, low_col)
      }

      reg_ls <- run_ls_series_nw(series_dt, panel, rhs)
      fmt_term_cell(reg_ls, tt, scale_est = scale_for(tt))
    })

    col_spread_mat <- matrix(
      col_spread,
      ncol = 1,
      dimnames = list(row_levels, col_spread_label)
    )

    base_mat <- cbind(base_mat, col_spread_mat)

    ## ------------------------------------------------------------
    ## Bottom row: spread across rows within each column
    ##
    ## Examples:
    ##   BM:          S1 - S5 within each B/M column
    ##   OP:          S1 - S5 within each OP column
    ##   INV:         S1 - S5 within each Inv column
    ##   CHI_SIC:     S1 - S5 within each chi column
    ##   BM_CHI_SIC:  B5 - B1 within each chi column
    ##   OP_CHI_SIC:  O5 - O1 within each chi column
    ##   INV_CHI_SIC: I1 - I5 within each chi column
    ## ------------------------------------------------------------

    row_spread <- sapply(col_levels, function(cc) {
      low_row_col <- paste0(spec$prefix, row_levels[1], cc)
      high_row_col <- paste0(spec$prefix, row_levels[5], cc)

      if (spec$row == "S") {
        ## Size spread is small-minus-big
        series_dt <- compute_ls_series_5x5(panel, low_row_col, high_row_col)
      } else if (spec$row == "I") {
        ## Investment spread is low-minus-high
        series_dt <- compute_ls_series_5x5(panel, low_row_col, high_row_col)
      } else {
        ## BM and OP use high-minus-low
        series_dt <- compute_ls_series_5x5(panel, high_row_col, low_row_col)
      }

      reg_ls <- run_ls_series_nw(series_dt, panel, rhs)
      fmt_term_cell(reg_ls, tt, scale_est = scale_for(tt))
    })

    ## Bottom-right corner intentionally left blank
    row_spread <- c(row_spread, "")

    row_spread_mat <- matrix(
      row_spread,
      nrow = 1,
      dimnames = list(row_spread_label, colnames(base_mat))
    )

    base_mat <- rbind(base_mat, row_spread_mat)

    result[[tt]] <- base_mat
  }

  result
}


### 8.9.3 Run and store all 25-portfolio loading tables #######################

## Main object with one entry per family-model combination
loading_tables_25 <- list()

## Same information, but nested by model.
## This makes it easier to print one model at a time.
loading_tables_25_by_model <- list()

for (mod in names(models)) {
  rhs <- models[[mod]]
  loading_tables_25_by_model[[mod]] <- list()

  for (fam in names(families_25)) {
    spec <- families_25[[fam]]

    result <- make_full_loading_table_25(
      panel      = factor_test_panel,
      fam        = fam,
      spec       = spec,
      model_name = mod,
      rhs        = rhs
    )

    key <- paste(fam, mod, sep = "_")

    loading_tables_25[[key]] <- result
    loading_tables_25_by_model[[mod]][[fam]] <- result
  }
}

saveRDS(loading_tables_25, "out_loading_tables_25port.rds")
saveRDS(loading_tables_25_by_model, "out_loading_tables_25port_by_model.rds")

### 8.9.4 Print 25-portfolio loading tables by model ##########################
print_loading_tables_25_model <- function(model_name) {
  if (!model_name %in% names(loading_tables_25_by_model)) {
    stop("Model not found: ", model_name)
  }

  cat("\n============================================================\n")
  cat("25-PORTFOLIO LOADING TABLES:", model_name, "\n")
  cat("============================================================\n")

  model_tables <- loading_tables_25_by_model[[model_name]]

  for (fam in names(model_tables)) {
    result <- model_tables[[fam]]
    if (is.null(result)) next

    cat("\n------------------------------------------------------------\n")
    cat("Family:", fam, "\n")
    cat("Model: ", model_name, "\n")
    cat("------------------------------------------------------------\n")

    for (tt in names(result)) {
      label <- if (tt == "(Intercept)") {
        "Alpha (%/month)"
      } else {
        tt
      }

      cat("\n", label, "\n", sep = "")
      print(result[[tt]])
    }
  }

  invisible(model_tables)
}


### 8.9.5 Print one model at a time ###########################################
print_loading_tables_25_model("CAPM")
print_loading_tables_25_model("FF3")
print_loading_tables_25_model("FF5")
print_loading_tables_25_model("FF5_SMQ")

### 8.9.6 Save full printed output to text files by model ######################
for (mod in names(models)) {
  file_name <- paste0("out_loading_tables_25port_", mod, ".txt")
  capture.output(
    print_loading_tables_25_model(mod),
    file = file_name
  )
  cat("Saved:", file_name, "\n")
}

### 8.10 32-Portfolio Families: Alphas and Loadings ###########################
## Structure mirrors the paper's Tables 8, 9, 10 for 2x4x4 sorts.
## For each family x model, one block per size group (Small / Big):
##   - Alpha panel, then each factor loading panel
##   - Rightmost column: spread across the column characteristic
##   - Bottom row: spread across the row characteristic
##
## Sign conventions:
##   - H-L for BM, OP, chi
##   - L-H for Inv
##
## For chi families:
##   - Rightmost column is H-L chi
##   - Bottom row is H-L or L-H characteristic spread
##
## For non-chi families:
##   - Rightmost column is the spread across the second characteristic
##   - Bottom row is the spread across the first characteristic
##
## The bottom-right corner is left blank because this would be a
## double-difference portfolio, which is not constructed here.


### 8.10.1 Long-short column mappings for 32-port chi families #################
## Kept because Section 8.11 uses this object.
## The new 8.10 helper below does not need these precomputed LS columns,
## because it constructs spreads directly from the 2x4x4 portfolio returns.

ls_chi_32_config <- list(
  OP_CHI_SIC = list(
    chi_pfx  = "LS_OP_CHI_SIC_LS_OP_CHI_SIC_",
    char_pfx = "LS_CHAR_OP_CHI_SIC_LS_OP_CHI_SIC_"
  ),
  BM_CHI_SIC = list(
    chi_pfx  = "LS_BM_CHI_SIC_LS_BM_CHI_SIC_",
    char_pfx = "LS_CHAR_BM_CHI_SIC_LS_BM_CHI_SIC_"
  ),
  INV_CHI_SIC = list(
    chi_pfx  = "LS_INV_CHI_SIC_LS_INV_CHI_SIC_",
    char_pfx = "LS_CHAR_INV_CHI_SIC_LS_INV_CHI_SIC_"
  )
)


### 8.10.2 Characteristic structure for 32-portfolio families #################
## c1 is the row characteristic.
## c2 is the column characteristic.

family_dims_32 <- list(
  BM_OP       = list(row_char = "BM", col_char = "OP"),
  BM_INV      = list(row_char = "BM", col_char = "INV"),
  OP_INV      = list(row_char = "OP", col_char = "INV"),
  OP_CHI_SIC  = list(row_char = "OP", col_char = "CHI"),
  BM_CHI_SIC  = list(row_char = "BM", col_char = "CHI"),
  INV_CHI_SIC = list(row_char = "INV", col_char = "CHI")
)

spread_label_32 <- function(char_name) {
  if (char_name == "CHI") {
    return("H-L chi")
  }
  if (char_name == "INV") {
    return("L-H char")
  }
  return("H-L char")
}

spread_long_short_order_32 <- function(char_name, low_col, high_col) {
  ## Returns list(long_col, short_col)
  ## BM, OP and chi: high-minus-low
  ## Inv: low-minus-high

  if (char_name == "INV") {
    return(list(long_col = low_col, short_col = high_col))
  }

  list(long_col = high_col, short_col = low_col)
}


### 8.10.3 Long-short helpers #################################################

## Build one long-short return series directly from two portfolio columns.
## Since both sides are portfolio returns, RF cancels in the long-short spread.
compute_ls_series_32 <- function(panel, long_col, short_col) {
  if (!long_col %in% names(panel) || !short_col %in% names(panel)) {
    return(NULL)
  }

  data.table(
    mdate = panel$mdate,
    ret   = panel[[long_col]] - panel[[short_col]]
  )
}

## Run NW regression on an already constructed long-short return series
run_ls_series_nw_32 <- function(series_dt, panel, rhs, lags = NW_LAGS) {
  if (is.null(series_dt)) {
    return(NULL)
  }

  dt <- merge(
    series_dt,
    panel[, c("mdate", rhs), with = FALSE],
    by = "mdate",
    all.x = TRUE
  )

  y <- dt$ret
  X <- dt[, ..rhs]

  run_ts_reg_nw(y, X, lags)
}

## Format any term from a long-short regression
fmt_term_cell_32 <- function(reg, term_name, scale_est = 1, d_est = 2, d_t = 2) {
  if (is.null(reg)) {
    return("")
  }

  b <- reg[term == term_name]
  if (nrow(b) == 0) {
    return("")
  }

  fmt_cell(
    est       = b$estimate,
    tstat     = b$tstat,
    stars     = b$stars,
    scale_est = scale_est,
    d_est     = d_est,
    d_t       = d_t
  )
}


### 8.10.4 Helper: build full loading table for one 2x4x4 family x model ######

## Returns:
##   list(
##     S = list(term = matrix),
##     B = list(term = matrix)
##   )

make_full_loading_table_32 <- function(panel, fam, spec, model_name, rhs) {
  reg <- run_regs_2x4x4(panel, spec$prefix, rhs)
  if (is.null(reg) || nrow(reg) == 0) {
    return(NULL)
  }

  if (!fam %in% names(family_dims_32)) {
    stop("Family dimensions not defined for: ", fam)
  }

  dim_cfg <- family_dims_32[[fam]]

  row_char <- dim_cfg$row_char
  col_char <- dim_cfg$col_char

  row_levels <- as.character(1:4)
  col_levels <- as.character(1:4)

  row_spread_label <- spread_label_32(row_char)
  col_spread_label <- spread_label_32(col_char)

  all_terms <- c("(Intercept)", rhs)

  ## Alpha is reported in monthly %, factor loadings are not scaled
  scale_for <- function(tt) ifelse(tt == "(Intercept)", 100, 1)

  result <- list(S = list(), B = list())

  for (sz in c("S", "B")) {
    for (tt in all_terms) {
      ## Base 4x4 coefficient matrix for this size group
      base_mat <- term_to_4x4(
        reg_dt     = reg,
        term_name  = tt,
        size_flag  = sz,
        scale_est  = scale_for(tt)
      )

      ## ------------------------------------------------------------
      ## Rightmost column: spread across columns within each row
      ##
      ## Examples:
      ##   BM_OP:       OP4 - OP1 within each BM row
      ##   BM_INV:      Inv1 - Inv4 within each BM row
      ##   OP_INV:      Inv1 - Inv4 within each OP row
      ##   OP_CHI_SIC:  Chi4 - Chi1 within each OP row
      ##   BM_CHI_SIC:  Chi4 - Chi1 within each B/M row
      ##   INV_CHI_SIC: Chi4 - Chi1 within each Inv row
      ## ------------------------------------------------------------

      col_spread <- sapply(row_levels, function(rr) {
        low_col_name <- paste0(spec$prefix, sz, rr, "1")
        high_col_name <- paste0(spec$prefix, sz, rr, "4")

        ls_order <- spread_long_short_order_32(
          char_name = col_char,
          low_col   = low_col_name,
          high_col  = high_col_name
        )

        series_dt <- compute_ls_series_32(
          panel     = panel,
          long_col  = ls_order$long_col,
          short_col = ls_order$short_col
        )

        reg_ls <- run_ls_series_nw_32(series_dt, panel, rhs)
        fmt_term_cell_32(reg_ls, tt, scale_est = scale_for(tt))
      })

      col_spread_mat <- matrix(
        col_spread,
        ncol = 1,
        dimnames = list(rownames(base_mat), col_spread_label)
      )

      base_mat <- cbind(base_mat, col_spread_mat)

      ## ------------------------------------------------------------
      ## Bottom row: spread across rows within each column
      ##
      ## Examples:
      ##   BM_OP:       BM4 - BM1 within each OP column
      ##   BM_INV:      BM4 - BM1 within each Inv column
      ##   OP_INV:      OP4 - OP1 within each Inv column
      ##   OP_CHI_SIC:  OP4 - OP1 within each chi column
      ##   BM_CHI_SIC:  BM4 - BM1 within each chi column
      ##   INV_CHI_SIC: Inv1 - Inv4 within each chi column
      ## ------------------------------------------------------------

      row_spread <- sapply(col_levels, function(cc) {
        low_row_col_name <- paste0(spec$prefix, sz, "1", cc)
        high_row_col_name <- paste0(spec$prefix, sz, "4", cc)

        ls_order <- spread_long_short_order_32(
          char_name = row_char,
          low_col   = low_row_col_name,
          high_col  = high_row_col_name
        )

        series_dt <- compute_ls_series_32(
          panel     = panel,
          long_col  = ls_order$long_col,
          short_col = ls_order$short_col
        )

        reg_ls <- run_ls_series_nw_32(series_dt, panel, rhs)
        fmt_term_cell_32(reg_ls, tt, scale_est = scale_for(tt))
      })

      ## Bottom-right corner intentionally left blank
      row_spread <- c(row_spread, "")

      row_spread_mat <- matrix(
        row_spread,
        nrow = 1,
        dimnames = list(row_spread_label, colnames(base_mat))
      )

      base_mat <- rbind(base_mat, row_spread_mat)

      result[[sz]][[tt]] <- base_mat
    }
  }

  result
}


### 8.10.5 Run and store all 32-portfolio loading tables ######################

## Main object with one entry per family-model combination
loading_tables_32 <- list()

## Same information, but nested by model.
## This makes it easier to print one model at a time.
loading_tables_32_by_model <- list()

for (mod in names(models)) {
  rhs <- models[[mod]]
  loading_tables_32_by_model[[mod]] <- list()

  for (fam in names(families_32)) {
    spec <- families_32[[fam]]

    result <- make_full_loading_table_32(
      panel      = factor_test_panel,
      fam        = fam,
      spec       = spec,
      model_name = mod,
      rhs        = rhs
    )

    key <- paste(fam, mod, sep = "_")

    loading_tables_32[[key]] <- result
    loading_tables_32_by_model[[mod]][[fam]] <- result
  }
}

saveRDS(loading_tables_32, "out_loading_tables_32port.rds")
saveRDS(loading_tables_32_by_model, "out_loading_tables_32port_by_model.rds")


### 8.10.6 Print 32-portfolio loading tables by model #########################

print_loading_tables_32_model <- function(model_name) {
  if (!model_name %in% names(loading_tables_32_by_model)) {
    stop("Model not found: ", model_name)
  }

  cat("\n============================================================\n")
  cat("32-PORTFOLIO LOADING TABLES:", model_name, "\n")
  cat("============================================================\n")

  model_tables <- loading_tables_32_by_model[[model_name]]

  for (fam in names(model_tables)) {
    result <- model_tables[[fam]]
    if (is.null(result)) next

    cat("\n------------------------------------------------------------\n")
    cat("Family:", fam, "\n")
    cat("Model: ", model_name, "\n")
    cat("------------------------------------------------------------\n")

    for (tt in names(result$S)) {
      label <- if (tt == "(Intercept)") {
        "Alpha (%/month)"
      } else {
        tt
      }

      cat("\n", label, "\n", sep = "")
      cat("Small\n")
      print(result$S[[tt]])
      cat("Big\n")
      print(result$B[[tt]])
    }
  }

  invisible(model_tables)
}


### 8.10.7 Print one model at a time ##########################################

## Run these manually when you want to inspect the console output.
print_loading_tables_32_model("CAPM")
print_loading_tables_32_model("FF3")
print_loading_tables_32_model("FF5")
print_loading_tables_32_model("FF5_SMQ")

### 8.10.8 Save full printed output to text files by model #####################
for (mod in names(models)) {
  file_name <- paste0("out_loading_tables_32port_", mod, ".txt")
  capture.output(
    print_loading_tables_32_model(mod),
    file = file_name
  )
  cat("Saved:", file_name, "\n")
}


### 8.11 Figures ##############################################################
## Uses FF5 regression results on all 25- and 32-portfolio test assets.
## Figures are printed directly to the RStudio Plots pane rather than saved.


### 8.11.1 Plot settings ######################################################
## Colour palettes, coefficient groups, family labels, and shared themes

## Five-group palette
size_cols_5 <- c(
  "S1" = "#08306B",
  "S2" = "#2171B5",
  "S3" = "#9ECAE1",
  "S4" = "#FC9272",
  "S5" = "#CB181D"
)

## Four-group palette
char_cols_4 <- c(
  "Q1" = "#08306B",
  "Q2" = "#6BAED6",
  "Q3" = "#FC9272",
  "Q4" = "#CB181D"
)

## Small/big palette
size_cols_2 <- c(
  "Small firms" = "#08306B",
  "Big firms"   = "#E1E3E5"
)

## Heatmap endpoints
heat_low <- "#08306B"
heat_high <- "#CB181D"

## Coefficient groups used in figures
coef_all <- c("Alpha", "MKT", "SMB", "HML", "RMW", "CMA")
coef_core <- c("Alpha", "HML", "RMW")

## 25-portfolio families
fams_25_plot <- c(
  "BM",
  "OP",
  "INV",
  "CHI_SIC",
  "BM_CHI_SIC",
  "OP_CHI_SIC",
  "INV_CHI_SIC"
)

## 32-portfolio families
fams_32_plot <- c(
  "BM_OP",
  "BM_INV",
  "OP_INV",
  "OP_CHI_SIC",
  "BM_CHI_SIC",
  "INV_CHI_SIC"
)

## 32-portfolio chi families used later in line and bar plots
chi_fams_32_plot <- c("OP_CHI_SIC", "BM_CHI_SIC", "INV_CHI_SIC")

## Shared theme for bar plots
pretty_bar_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey35"),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(colour = "grey20"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 10),
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F3F4F6", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

## Clean family labels for 25-portfolio figures
pretty_family_label_25 <- function(x) {
  fcase(
    x == "BM",          "Size-B/M",
    x == "OP",          "Size-OP",
    x == "INV",         "Size-Inv",
    x == "CHI_SIC",     "Size-chi",
    x == "BM_CHI_SIC",  "B/M-chi",
    x == "OP_CHI_SIC",  "OP-chi",
    x == "INV_CHI_SIC", "Inv-chi",
    default = x
  )
}

## Clean family labels for 32-portfolio figures
pretty_family_label_32 <- function(x) {
  fcase(
    x == "BM_OP",       "Size-B/M-OP",
    x == "BM_INV",      "Size-B/M-Inv",
    x == "OP_INV",      "Size-OP-Inv",
    x == "OP_CHI_SIC",  "Size-OP-chi",
    x == "BM_CHI_SIC",  "Size-B/M-chi",
    x == "INV_CHI_SIC", "Size-Inv-chi",
    default = x
  )
}

## Characteristic axis labels for 32-portfolio figures
row_axis_lbl_32 <- function(fam) {
  fcase(
    fam == "BM_OP",       "B/M quartile",
    fam == "BM_INV",      "B/M quartile",
    fam == "OP_INV",      "OP quartile",
    fam == "OP_CHI_SIC",  "OP quartile",
    fam == "BM_CHI_SIC",  "B/M quartile",
    fam == "INV_CHI_SIC", "Investment quartile",
    default = "First characteristic quartile"
  )
}

col_axis_lbl_32 <- function(fam) {
  fcase(
    fam == "BM_OP",       "OP quartile",
    fam == "BM_INV",      "Investment quartile",
    fam == "OP_INV",      "Investment quartile",
    fam == "OP_CHI_SIC",  "Chi quartile",
    fam == "BM_CHI_SIC",  "Chi quartile",
    fam == "INV_CHI_SIC", "Chi quartile",
    default = "Second characteristic quartile"
  )
}

## Characteristic labels used later in line plots
char_axis_lbl <- function(fam) {
  fcase(
    fam == "OP_CHI_SIC",  "OP quartile",
    fam == "BM_CHI_SIC",  "B/M quartile",
    fam == "INV_CHI_SIC", "Investment quartile",
    default = "Characteristic quartile"
  )
}


### 8.11.2 Portfolio-label helpers ############################################
## Extract family, size, and bucket labels from portfolio column names

label_portfolio <- function(port) {
  ## 25-portfolio patterns: Size x characteristic
  if (grepl("^BM_S[1-5]V[1-5]$", port)) {
    return(list(
      family   = "BM",
      size     = sub("^BM_(S\\d)V\\d$", "\\1", port),
      bucket_1 = sub("^BM_S\\d(V\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  if (grepl("^OP_S[1-5]P[1-5]$", port)) {
    return(list(
      family   = "OP",
      size     = sub("^OP_(S\\d)P\\d$", "\\1", port),
      bucket_1 = sub("^OP_S\\d(P\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  if (grepl("^INV_S[1-5]I[1-5]$", port)) {
    return(list(
      family   = "INV",
      size     = sub("^INV_(S\\d)I\\d$", "\\1", port),
      bucket_1 = sub("^INV_S\\d(I\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  if (grepl("^CHI_SIC_S[1-5]C[1-5]$", port)) {
    return(list(
      family   = "CHI_SIC",
      size     = sub("^CHI_SIC_(S\\d)C\\d$", "\\1", port),
      bucket_1 = sub("^CHI_SIC_S\\d(C\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  ## 25-portfolio patterns: characteristic x chi
  if (grepl("^BM_CHI_SIC_B[1-5]C[1-5]$", port)) {
    return(list(
      family   = "BM_CHI_SIC",
      size     = sub("^BM_CHI_SIC_(B\\d)C\\d$", "\\1", port),
      bucket_1 = sub("^BM_CHI_SIC_B\\d(C\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  if (grepl("^OP_CHI_SIC_O[1-5]C[1-5]$", port)) {
    return(list(
      family   = "OP_CHI_SIC",
      size     = sub("^OP_CHI_SIC_(O\\d)C\\d$", "\\1", port),
      bucket_1 = sub("^OP_CHI_SIC_O\\d(C\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  if (grepl("^INV_CHI_SIC_I[1-5]C[1-5]$", port)) {
    return(list(
      family   = "INV_CHI_SIC",
      size     = sub("^INV_CHI_SIC_(I\\d)C\\d$", "\\1", port),
      bucket_1 = sub("^INV_CHI_SIC_I\\d(C\\d)$", "\\1", port),
      bucket_2 = NA
    ))
  }

  ## 32-portfolio patterns: Size x B/M x OP
  if (grepl("^BM_OP_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "BM_OP",
      size     = sub("^BM_OP_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^BM_OP_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^BM_OP_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  if (grepl("^BM_OP_2X4X4_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "BM_OP",
      size     = sub("^BM_OP_2X4X4_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^BM_OP_2X4X4_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^BM_OP_2X4X4_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## 32-portfolio patterns: Size x B/M x Investment
  if (grepl("^BM_INV_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "BM_INV",
      size     = sub("^BM_INV_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^BM_INV_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^BM_INV_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  if (grepl("^BM_INV_2X4X4_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "BM_INV",
      size     = sub("^BM_INV_2X4X4_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^BM_INV_2X4X4_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^BM_INV_2X4X4_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## 32-portfolio patterns: Size x OP x Investment
  if (grepl("^OP_INV_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "OP_INV",
      size     = sub("^OP_INV_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^OP_INV_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^OP_INV_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  if (grepl("^OP_INV_2X4X4_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "OP_INV",
      size     = sub("^OP_INV_2X4X4_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^OP_INV_2X4X4_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^OP_INV_2X4X4_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## 32-portfolio patterns: Size x OP x chi
  if (grepl("^OP_CHI_SIC_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "OP_CHI_SIC",
      size     = sub("^OP_CHI_SIC_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^OP_CHI_SIC_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^OP_CHI_SIC_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## 32-portfolio patterns: Size x B/M x chi
  if (grepl("^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "BM_CHI_SIC",
      size     = sub("^BM_CHI_SIC_2X4X4_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^BM_CHI_SIC_2X4X4_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^BM_CHI_SIC_2X4X4_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## 32-portfolio patterns: Size x Investment x chi
  if (grepl("^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", port)) {
    return(list(
      family   = "INV_CHI_SIC",
      size     = sub("^INV_CHI_SIC_2X4X4_([SB])\\d\\d$", "\\1", port),
      bucket_1 = paste0("Q", sub("^INV_CHI_SIC_2X4X4_[SB](\\d)\\d$", "\\1", port)),
      bucket_2 = paste0("Q", sub("^INV_CHI_SIC_2X4X4_[SB]\\d(\\d)$", "\\1", port))
    ))
  }

  ## Default return if no pattern is matched
  list(family = NA, size = NA, bucket_1 = NA, bucket_2 = NA)
}


### 8.11.3 Rebuild FF5 regression results #####################################
## Re-estimate FF5 regressions for all plotted test assets

rhs_ff5 <- models[["FF5"]]

run_ff5_all <- function(panel, asset_cols, rhs) {
  ## Extract FF5 regressors
  X <- panel[, ..rhs]

  res <- lapply(asset_cols, function(col) {
    ## Portfolio excess return
    y <- panel[[col]] - panel[["RF"]]

    ## Time-series regression with Newey-West standard errors
    est <- run_ts_reg_nw(y, X)

    if (is.null(est)) {
      return(NULL)
    }

    ## Add portfolio metadata
    lab <- label_portfolio(col)

    est[, `:=`(
      portfolio = col,
      family    = lab$family,
      size      = lab$size,
      bucket_1  = lab$bucket_1,
      bucket_2  = lab$bucket_2
    )]

    est
  })

  ## Stack regression outputs
  rbindlist(res, fill = TRUE)
}

regs_ff5_25 <- run_ff5_all(factor_test_panel, all_25_cols, rhs_ff5)
regs_ff5_32 <- run_ff5_all(factor_test_panel, all_32_cols, rhs_ff5)


### 8.11.4 Prepare plotting data ##############################################
## Clean coefficient names, scale alphas, and add significance stars

make_plot_dt <- function(reg_dt) {
  out <- copy(reg_dt)

  ## Rename intercept and scale alpha to monthly percent
  out[, coefficient := fifelse(term == "(Intercept)", "Alpha", term)]
  out[, value := fifelse(term == "(Intercept)", 100 * estimate, estimate)]

  ## Add significance stars
  out[, sig := star_p(pval)]

  out[]
}

plot_dt_25 <- make_plot_dt(regs_ff5_25)
plot_dt_32 <- make_plot_dt(regs_ff5_32)


### 8.11.5 Plot helper functions ##############################################
## Shared helpers for heatmaps, line plots, and bar plots

## Diagnostic output: check which families were detected
cat("\nDetected 25-portfolio families:\n")
print(plot_dt_25[, .N, by = family][order(family)])

cat("\nDetected 32-portfolio families:\n")
print(plot_dt_32[, .N, by = family][order(family)])

## Normalise fill values within each coefficient panel
add_heatmap_fill <- function(dt) {
  out <- copy(dt)

  out[, max_abs := max(abs(value), na.rm = TRUE), by = coefficient]
  out[is.na(max_abs) | max_abs == 0, max_abs := 1]
  out[, fill_value := value / max_abs]

  out[]
}

## Row labels for 25-portfolio heatmaps
row_axis_lbl_25 <- function(fam) {
  fcase(
    fam == "BM",          "Size quintile",
    fam == "OP",          "Size quintile",
    fam == "INV",         "Size quintile",
    fam == "CHI_SIC",     "Size quintile",
    fam == "BM_CHI_SIC",  "B/M quintile",
    fam == "OP_CHI_SIC",  "OP quintile",
    fam == "INV_CHI_SIC", "Investment quintile",
    default = "Row portfolio"
  )
}

## Column labels for 25-portfolio heatmaps
col_axis_lbl_25 <- function(fam) {
  fcase(
    fam == "BM",          "B/M quintile",
    fam == "OP",          "OP quintile",
    fam == "INV",         "Investment quintile",
    fam == "CHI_SIC",     "Chi quintile",
    fam == "BM_CHI_SIC",  "Chi quintile",
    fam == "OP_CHI_SIC",  "Chi quintile",
    fam == "INV_CHI_SIC", "Chi quintile",
    default = "Column portfolio"
  )
}

## Row levels for 25-portfolio heatmaps
row_levels_25 <- function(fam) {
  if (fam %in% c("BM", "OP", "INV", "CHI_SIC")) {
    return(paste0("S", 1:5))
  }

  if (fam == "BM_CHI_SIC") {
    return(paste0("B", 1:5))
  }

  if (fam == "OP_CHI_SIC") {
    return(paste0("O", 1:5))
  }

  if (fam == "INV_CHI_SIC") {
    return(paste0("I", 1:5))
  }

  unique(plot_dt_25[family == fam, size])
}

## Column levels for 25-portfolio heatmaps
col_levels_25 <- function(fam) {
  if (fam == "BM") {
    return(paste0("V", 1:5))
  }

  if (fam == "OP") {
    return(paste0("P", 1:5))
  }

  if (fam == "INV") {
    return(paste0("I", 1:5))
  }

  if (fam %in% c("CHI_SIC", "BM_CHI_SIC", "OP_CHI_SIC", "INV_CHI_SIC")) {
    return(paste0("C", 1:5))
  }

  unique(plot_dt_25[family == fam, bucket_1])
}

## Characteristic labels for 32-portfolio line plots
char_labels_32 <- function(fam) {
  switch(fam,
    OP_CHI_SIC = c(
      "Q1" = "Q1 Low OP",
      "Q2" = "Q2",
      "Q3" = "Q3",
      "Q4" = "Q4 High OP"
    ),
    BM_CHI_SIC = c(
      "Q1" = "Q1 Low B/M",
      "Q2" = "Q2",
      "Q3" = "Q3",
      "Q4" = "Q4 High B/M"
    ),
    INV_CHI_SIC = c(
      "Q1" = "Q1 Low Inv",
      "Q2" = "Q2",
      "Q3" = "Q3",
      "Q4" = "Q4 High Inv"
    ),
    c("Q1" = "Q1", "Q2" = "Q2", "Q3" = "Q3", "Q4" = "Q4")
  )
}

## Coefficients shown for each 32-portfolio family in selected line plots
coefs_for_fam <- list(
  OP_CHI_SIC  = c("Alpha", "HML", "RMW"),
  BM_CHI_SIC  = c("Alpha", "HML", "RMW", "CMA"),
  INV_CHI_SIC = c("Alpha", "HML", "RMW")
)


### 8.11.6 Heatmaps ###########################################################
## FF5 alpha and factor-loading heatmaps.
## Cell values are actual coefficients; colours are scaled within each coefficient panel.


### 8.11.6.1 Heatmap helper for 25-portfolio families #########################
## Rows and columns refer to the two sorting dimensions in each 5 x 5 family.

plot_heatmap_25 <- function(fam) {
  dt <- copy(plot_dt_25[family == fam & coefficient %in% coef_all])

  if (nrow(dt) == 0) {
    stop(paste0("No observations found for 25-portfolio family: ", fam))
  }

  dt[, coefficient := factor(coefficient, levels = coef_all)]
  dt[, Row := factor(size, levels = row_levels_25(fam))]
  dt[, Col := factor(bucket_1, levels = col_levels_25(fam))]

  dt <- add_heatmap_fill(dt)

  ggplot(dt, aes(x = Col, y = Row, fill = fill_value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f%s", value, sig)), size = 2.6) +
    scale_fill_gradient2(
      low      = heat_low,
      mid      = "white",
      high     = heat_high,
      midpoint = 0,
      limits   = c(-1, 1)
    ) +
    facet_wrap(~coefficient, ncol = 3) +
    labs(
      title    = paste0(pretty_family_label_25(fam), ": FF5 alpha and factor loadings"),
      subtitle = "5 x 5 portfolios. Cell values are actual coefficients; colours scaled within each panel",
      x        = col_axis_lbl_25(fam),
      y        = row_axis_lbl_25(fam),
      fill     = "Scaled value"
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
}

### 8.11.6.2 Heatmaps for 25-portfolio families ###############################
## Each figure below refers to one 5 x 5 portfolio family.

## Size x B/M portfolios
cat("\nPrinting 25-portfolio heatmap: Size x B/M\n")
p_heat_25_bm <- plot_heatmap_25("BM")
print(p_heat_25_bm)

## Size x OP portfolios
cat("\nPrinting 25-portfolio heatmap: Size x OP\n")
p_heat_25_op <- plot_heatmap_25("OP")
print(p_heat_25_op)

## Size x Investment portfolios
cat("\nPrinting 25-portfolio heatmap: Size x Investment\n")
p_heat_25_inv <- plot_heatmap_25("INV")
print(p_heat_25_inv)

## Size x chi portfolios
cat("\nPrinting 25-portfolio heatmap: Size x chi\n")
p_heat_25_chi <- plot_heatmap_25("CHI_SIC")
print(p_heat_25_chi)

## B/M x chi portfolios
cat("\nPrinting 25-portfolio heatmap: B/M x chi\n")
p_heat_25_bm_chi <- plot_heatmap_25("BM_CHI_SIC")
print(p_heat_25_bm_chi)

## OP x chi portfolios
cat("\nPrinting 25-portfolio heatmap: OP x chi\n")
p_heat_25_op_chi <- plot_heatmap_25("OP_CHI_SIC")
print(p_heat_25_op_chi)

## Investment x chi portfolios
cat("\nPrinting 25-portfolio heatmap: Investment x chi\n")
p_heat_25_inv_chi <- plot_heatmap_25("INV_CHI_SIC")
print(p_heat_25_inv_chi)


### 8.11.6.3 Heatmap helper for 32-portfolio families #########################
## Rows refer to the first non-size characteristic.
## Columns refer to the second non-size characteristic.
## Panels split estimates by small and big firms.

plot_heatmap_32 <- function(fam) {
  dt <- copy(plot_dt_32[family == fam & coefficient %in% coef_all])

  if (nrow(dt) == 0) {
    stop(paste0("No observations found for 32-portfolio family: ", fam))
  }

  dt[, coefficient := factor(coefficient, levels = coef_all)]
  dt[, SizeGroup := fifelse(size == "S", "Small firms", "Big firms")]
  dt[, SizeGroup := factor(SizeGroup, levels = c("Small firms", "Big firms"))]
  dt[, Row := factor(bucket_1, levels = paste0("Q", 1:4))]
  dt[, Col := factor(bucket_2, levels = paste0("Q", 1:4))]

  dt <- add_heatmap_fill(dt)

  ggplot(dt, aes(x = Col, y = Row, fill = fill_value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f%s", value, sig)), size = 2.3) +
    scale_fill_gradient2(
      low      = heat_low,
      mid      = "white",
      high     = heat_high,
      midpoint = 0,
      limits   = c(-1, 1)
    ) +
    facet_grid(coefficient ~ SizeGroup) +
    labs(
      title    = paste0(pretty_family_label_32(fam), ": FF5 alpha and factor loadings"),
      subtitle = "2 x 4 x 4 portfolios. Cell values are actual coefficients; colours scaled within each panel",
      x        = col_axis_lbl_32(fam),
      y        = row_axis_lbl_32(fam),
      fill     = "Scaled value"
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "bottom"
    )
}


### 8.11.6.4 Heatmaps for 32-portfolio families ###############################
## Each figure below refers to one 2 x 4 x 4 portfolio family.

## Size x B/M x OP portfolios
cat("\nPrinting 32-portfolio heatmap: Size x B/M x OP\n")
p_heat_32_bm_op <- plot_heatmap_32("BM_OP")
print(p_heat_32_bm_op)

## Size x B/M x Investment portfolios
cat("\nPrinting 32-portfolio heatmap: Size x B/M x Investment\n")
p_heat_32_bm_inv <- plot_heatmap_32("BM_INV")
print(p_heat_32_bm_inv)

## Size x OP x Investment portfolios
cat("\nPrinting 32-portfolio heatmap: Size x OP x Investment\n")
p_heat_32_op_inv <- plot_heatmap_32("OP_INV")
print(p_heat_32_op_inv)

## Size x OP x chi portfolios
cat("\nPrinting 32-portfolio heatmap: Size x OP x chi\n")
p_heat_32_op_chi <- plot_heatmap_32("OP_CHI_SIC")
print(p_heat_32_op_chi)

## Size x B/M x chi portfolios
cat("\nPrinting 32-portfolio heatmap: Size x B/M x chi\n")
p_heat_32_bm_chi <- plot_heatmap_32("BM_CHI_SIC")
print(p_heat_32_bm_chi)

## Size x Investment x chi portfolios
cat("\nPrinting 32-portfolio heatmap: Size x Investment x chi\n")
p_heat_32_inv_chi <- plot_heatmap_32("INV_CHI_SIC")
print(p_heat_32_inv_chi)


### 8.11.7 Line plots #########################################################
## Alpha, HML, and RMW profiles across chi portfolios.
## Figures are printed separately for each relevant portfolio family.


### 8.11.7.1 Line plot helper for 25-portfolio Size x chi ######################
## Shows coefficient profiles across chi quintiles, separately by size quintile.

plot_line_25_chi <- function() {
  dt <- copy(plot_dt_25[
    family == "CHI_SIC" & coefficient %in% coef_core
  ])

  if (nrow(dt) == 0) {
    stop("No observations found for 25-portfolio family: CHI_SIC")
  }

  ## Format variables
  dt[, coefficient := factor(coefficient, levels = coef_core)]
  dt[, Size := factor(size, levels = paste0("S", 1:5))]
  dt[, Chi := factor(bucket_1, levels = paste0("C", 1:5))]

  ## Plot
  ggplot(dt, aes(x = Chi, y = value, group = Size, color = Size)) +
    geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.4
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_wrap(
      ~coefficient,
      ncol   = 1,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = size_cols_5,
      labels = c(
        "S1" = "S1 Smallest",
        "S2" = "S2",
        "S3" = "S3",
        "S4" = "S4",
        "S5" = "S5 Largest"
      )
    ) +
    labs(
      title    = "Size-chi: FF5 alpha, HML, and RMW profiles",
      subtitle = "Profiles across chi quintiles, separately by size quintile",
      x        = "Chi quintile",
      y        = "Coefficient value",
      color    = "Size quintile"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 10, colour = "grey35"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}


### 8.11.7.2 Line plot helper for 32-portfolio chi families ####################
## Shows coefficient profiles across chi quartiles.
## Lines represent the non-chi characteristic quartiles.
## Panels separate coefficients and size groups.

plot_line_32_chi <- function(fam) {
  if (!fam %in% names(coefs_for_fam)) {
    stop("No coefficient selection defined for family: ", fam)
  }

  dt <- copy(plot_dt_32[
    family == fam & coefficient %in% coefs_for_fam[[fam]]
  ])

  if (nrow(dt) == 0) {
    stop(paste0("No observations found for 32-portfolio family: ", fam))
  }

  ## Format variables
  dt[, coefficient := factor(coefficient, levels = coefs_for_fam[[fam]])]
  dt[, SizeGroup := fifelse(size == "S", "Small firms", "Big firms")]
  dt[, SizeGroup := factor(SizeGroup, levels = c("Small firms", "Big firms"))]
  dt[, Chi := factor(bucket_2, levels = paste0("Q", 1:4))]
  dt[, Char := factor(bucket_1, levels = paste0("Q", 1:4))]

  ## Characteristic labels depend on the portfolio family
  char_labs <- char_labels_32(fam)

  ## Plot
  ggplot(dt, aes(x = Chi, y = value, group = Char, color = Char)) +
    geom_hline(
      yintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.4
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(
      coefficient ~ SizeGroup,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = char_cols_4,
      labels = char_labs
    ) +
    labs(
      title    = paste0(pretty_family_label_32(fam), ": FF5 coefficient profiles across chi"),
      subtitle = "Lines show profiles across chi quartiles within each characteristic quartile",
      x        = "Chi quartile",
      y        = "Coefficient value",
      color    = char_axis_lbl(fam)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 10, colour = "grey35"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}


### 8.11.7.3 Line plot for 25-portfolio Size x chi ############################
## Figure refers to the 5 x 5 Size x chi portfolio family.

cat("\nPrinting 25-portfolio line plot: Size x chi\n")
p_line_25_chi <- plot_line_25_chi()
print(p_line_25_chi)


### 8.11.7.4 Line plots for 32-portfolio chi families #########################
## Each figure below refers to one 2 x 4 x 4 portfolio family.

## Size x OP x chi portfolios
cat("\nPrinting 32-portfolio line plot: Size x OP x chi\n")
p_line_32_op_chi <- plot_line_32_chi("OP_CHI_SIC")
print(p_line_32_op_chi)

## Size x B/M x chi portfolios
cat("\nPrinting 32-portfolio line plot: Size x B/M x chi\n")
p_line_32_bm_chi <- plot_line_32_chi("BM_CHI_SIC")
print(p_line_32_bm_chi)

## Size x Investment x chi portfolios
cat("\nPrinting 32-portfolio line plot: Size x Investment x chi\n")
p_line_32_inv_chi <- plot_line_32_chi("INV_CHI_SIC")
print(p_line_32_inv_chi)

### 8.11.8 Bar plots ##########################################################
## Small-versus-big averages and high-minus-low chi spreads

### 8.11.8.1 Average coefficients by chi quartile #############################
## Average across the non-chi characteristic within each 32-portfolio chi family.

summary_by_chi <- copy(plot_dt_32[
  family %in% chi_fams_32_plot & coefficient %in% coef_core,
  .(value = mean(value, na.rm = TRUE)),
  by = .(family, coefficient, size, bucket_2)
])

## Add clean family labels for plotting
summary_by_chi[, Family := factor(
  pretty_family_label_32(family),
  levels = sapply(chi_fams_32_plot, pretty_family_label_32)
)]

## Format grouping variables for the plot
summary_by_chi[, coefficient := factor(coefficient, levels = coef_core)]
summary_by_chi[, SizeGroup := fifelse(size == "S", "Small firms", "Big firms")]
summary_by_chi[, SizeGroup := factor(SizeGroup, levels = c("Small firms", "Big firms"))]
summary_by_chi[, Chi := factor(bucket_2, levels = paste0("Q", 1:4))]

## Plot average coefficients by chi quartile and size group
p_bars_chi <- ggplot(summary_by_chi, aes(x = Chi, y = value, fill = SizeGroup)) +
  geom_col(
    position = position_dodge(0.75),
    width    = 0.62,
    colour   = NA,
    alpha    = 0.95
  ) +
  facet_grid(coefficient ~ Family, scales = "free_y") +
  scale_fill_manual(values = size_cols_2) +
  labs(
    title    = "Small versus Big Firms: Average Alpha, HML, and RMW by Chi Quartile",
    subtitle = "Averages across OP, B/M, or investment quartiles within each family",
    x        = "Chi quartile",
    y        = "Average coefficient value",
    fill     = "Size group"
  ) +
  pretty_bar_theme

print(p_bars_chi)

ggsave(
  filename = file.path(FIGURES_DIR, "figure_small_big_avg_coefficients_by_chi.pdf"),
  plot     = p_bars_chi,
  width    = 12,
  height   = 12
)

normalizePath(file.path(FIGURES_DIR, "figure_small_big_avg_coefficients_by_chi.pdf"))

### 8.11.8.2 High-minus-low chi spreads by characteristic quartile #############
## Spread is computed as the highest chi quartile minus the lowest chi quartile.

spread_input <- copy(plot_dt_32[
  family %in% chi_fams_32_plot & coefficient %in% coef_core,
  .(family, coefficient, size, bucket_1, bucket_2, value)
])

## Reshape chi quartiles to columns
spread_wide <- dcast(
  spread_input,
  family + coefficient + size + bucket_1 ~ bucket_2,
  value.var = "value"
)

## Compute high-minus-low chi spread
spread_wide[, chi_spread := Q4 - Q1]

## Add clean family labels for plotting
spread_wide[, Family := factor(
  pretty_family_label_32(family),
  levels = sapply(chi_fams_32_plot, pretty_family_label_32)
)]

## Format grouping variables for the plot
spread_wide[, coefficient := factor(coefficient, levels = coef_core)]
spread_wide[, SizeGroup := fifelse(size == "S", "Small firms", "Big firms")]
spread_wide[, SizeGroup := factor(SizeGroup, levels = c("Small firms", "Big firms"))]
spread_wide[, Characteristic := factor(bucket_1, levels = paste0("Q", 1:4))]

## Plot high-minus-low chi spreads by characteristic quartile and size group
p_bars_spread <- ggplot(
  spread_wide,
  aes(x = Characteristic, y = chi_spread, fill = SizeGroup)
) +
  geom_col(
    position = position_dodge(0.75),
    width    = 0.62,
    colour   = NA,
    alpha    = 0.95
  ) +
  facet_grid(coefficient ~ Family, scales = "free_y") +
  scale_fill_manual(values = size_cols_2) +
  labs(
    title    = "Small versus Big Firms: High-minus-Low Chi Spreads",
    subtitle = "Spread = highest chi quartile minus lowest chi quartile",
    x        = "OP, B/M, or Inv quartile",
    y        = "High-chi minus low-chi coefficient",
    fill     = "Size group"
  ) +
  pretty_bar_theme

print(p_bars_spread)


### 8.11.8.3 Average high-minus-low chi spreads ###############################
## Average the chi spread across characteristic quartiles.
avg_spread <- spread_wide[
  ,
  .(avg_chi_spread = mean(chi_spread, na.rm = TRUE)),
  by = .(Family, coefficient, SizeGroup)
]

## Format coefficient order for plotting
avg_spread[, coefficient := factor(coefficient, levels = coef_core)]

## Plot average high-minus-low chi spreads by family and size group
p_bars_avg <- ggplot(
  avg_spread,
  aes(x = Family, y = avg_chi_spread, fill = SizeGroup)
) +
  geom_col(
    position = position_dodge(0.75),
    width    = 0.62,
    colour   = NA,
    alpha    = 0.95
  ) +
  facet_wrap(~coefficient, scales = "free_y") +
  scale_fill_manual(values = size_cols_2) +
  labs(
    title    = "Small versus Big Firms: Average High-minus-Low Chi Spreads",
    subtitle = "Average across OP, B/M, or investment quartiles within each family",
    x        = "",
    y        = "Average high-chi minus low-chi coefficient",
    fill     = "Size group"
  ) +
  pretty_bar_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

print(p_bars_avg)




### Pairwise Spanning Tests: HML, CMA, and RMW ###############################

factor_pairs <- list(
  HML_CMA = c("HML", "CMA"),
  HML_RMW = c("HML", "RMW"),
  CMA_RMW = c("CMA", "RMW")
)

pairwise_spanning_results <- lapply(
  factor_pairs,
  function(pair) {
    make_spanning_table(
      panel      = factor_test_panel,
      factor_vec = pair
    )
  }
)

cat("\n============================================================\n")
cat("PAIRWISE SPANNING REGRESSIONS: HML, CMA, AND RMW\n")
cat("============================================================\n")

for (pair_name in names(pairwise_spanning_results)) {
  cat("\n---", pair_name, "---\n")
  print(pairwise_spanning_results[[pair_name]])
}

saveRDS(
  pairwise_spanning_results,
  "out_pairwise_spanning_HML_CMA_RMW.rds"
)
