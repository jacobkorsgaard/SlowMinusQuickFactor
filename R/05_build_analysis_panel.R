## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026

### 5. Time Series Regressions ################################################
### 5.0 Combined Setup ########################################################

## Clear Everything
cat("\014")
rm(list = ls())
graphics.off()

## Load Libraries
library(data.table)
library(lubridate)
library(zoo)
library(sandwich)
library(lmtest)
library(xtable)

## Define Window
START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")

## Newey-West SE Lags
NW_LAGS <- 12

## Project paths
source(file.path("R", "00_config.R"))

### Helper Functions ##########################################################
ensure_mdate <- function(dt) {
  setDT(dt)
  dt[, mdate := as.Date(mdate)]
  dt
}

load_and_prefix <- function(file, prefix) {
  dt <- readRDS(file)
  setDT(dt)
  dt[, mdate := as.Date(mdate)]
  cols <- setdiff(names(dt), "mdate")
  setnames(dt, cols, paste0(prefix, cols))
  dt
}

safe_share <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

### 5.0.1 Load Factor Data ####################################################

## FF5 factors
ff5 <- readRDS("ff5_factors_internal.rds")
setDT(ff5)
ff5[, mdate := as.Date(mdate)]

## SMQ factor
SMQ <- readRDS("SMQ_factors.rds")
setDT(SMQ)
SMQ[, mdate := as.Date(mdate)]

## Risk-free rate source
rf_src <- readRDS("crsp_compustat_with_adj_cost.rds")
setDT(rf_src)
rf_src[, mdate := as.Date(mdate)]

rf_monthly <-
  unique(
    rf_src[
      !is.na(RF),
      .(mdate, RF)
    ]
  )

setkey(rf_monthly, mdate)

### 5.0.2 Load 1D Test Assets #################################################
size_1d <- load_and_prefix("my_SIZE_1D_wide.rds", "SIZE_")
bm1d_1d <- load_and_prefix("my_BM1D_1D_wide.rds", "BM1D_")
op1d_1d <- load_and_prefix("my_OP1D_1D_wide.rds", "OP1D_")
inv1d_1d <- load_and_prefix("my_INV1D_1D_wide.rds", "INV1D_")
chi1d_1d <- load_and_prefix("my_CHI1D_1D_wide.rds", "CHI1D_")

### 5.0.3 Load 5x5 Test Assets ################################################
bm_5x5 <- load_and_prefix("my_Size_BM_5x5_wide.rds", "BM_")
op_5x5 <- load_and_prefix("my_Size_OP_5x5_wide.rds", "OP_")
inv_5x5 <- load_and_prefix("my_Size_INV_5x5_wide.rds", "INV_")
chi_firm_5x5 <- load_and_prefix("my_Size_chi_firm_5x5_wide.rds", "CHI_FIRM_")
chi_firm_exp_5x5 <- load_and_prefix("my_Size_chi_firm_exp_5x5_wide.rds", "CHI_FIRM_EXP_")
chi_sic_5x5 <- load_and_prefix("my_Size_chi_sic_5x5_wide.rds", "CHI_SIC_")
chi_sic_exp_5x5 <- load_and_prefix("my_Size_chi_sic_exp_5x5_wide.rds", "CHI_SIC_EXP_")
chi_sic2_5x5 <- load_and_prefix("my_Size_chi_sic2_5x5_wide.rds", "CHI_SIC2_")
chi_sic2_exp_5x5 <- load_and_prefix("my_Size_chi_sic2_exp_5x5_wide.rds", "CHI_SIC2_EXP_")

## New 5x5 chi interaction portfolios
bm_chi_sic_5x5 <- load_and_prefix("my_BM_CHI_SIC_5x5_wide.rds", "BM_CHI_SIC_")
op_chi_sic_5x5 <- load_and_prefix("my_OP_CHI_SIC_5x5_wide.rds", "OP_CHI_SIC_")
inv_chi_sic_5x5 <- load_and_prefix("my_INV_CHI_SIC_5x5_wide.rds", "INV_CHI_SIC_")

### 5.0.4 Load 2x4x4 Test Assets ##############################################
bm_op_2x4x4 <- load_and_prefix("my_Size_BM_OP_2x4x4_wide.rds", "BM_OP_")
bm_inv_2x4x4 <- load_and_prefix("my_Size_BM_INV_2x4x4_wide.rds", "BM_INV_")
op_inv_2x4x4 <- load_and_prefix("my_Size_OP_INV_2x4x4_wide.rds", "OP_INV_")

op_chi_firm_2x4x4 <- load_and_prefix("my_Size_OP_CHI_FIRM_2x4x4_wide.rds", "OP_CHI_FIRM_")
op_chi_sic_2x4x4 <- load_and_prefix("my_Size_OP_CHI_SIC_2x4x4_wide.rds", "OP_CHI_SIC_")
op_chi_sic2_2x4x4 <- load_and_prefix("my_Size_OP_CHI_SIC2_2x4x4_wide.rds", "OP_CHI_SIC2_")

bm_chi_firm_2x4x4 <- load_and_prefix("my_Size_BM_CHI_FIRM_2x4x4_wide.rds", "BM_CHI_FIRM_")
bm_chi_sic_2x4x4 <- load_and_prefix("my_Size_BM_CHI_SIC_2x4x4_wide.rds", "BM_CHI_SIC_2X4X4_")
bm_chi_sic2_2x4x4 <- load_and_prefix("my_Size_BM_CHI_SIC2_2x4x4_wide.rds", "BM_CHI_SIC2_")

inv_chi_firm_2x4x4 <- load_and_prefix("my_Size_INV_CHI_FIRM_2x4x4_wide.rds", "INV_CHI_FIRM_")
inv_chi_sic_2x4x4 <- load_and_prefix("my_Size_INV_CHI_SIC_2x4x4_wide.rds", "INV_CHI_SIC_2X4X4_")
inv_chi_sic2_2x4x4 <- load_and_prefix("my_Size_INV_CHI_SIC2_2x4x4_wide.rds", "INV_CHI_SIC2_")

### 5.0.5 Load Chi Decile (1x10) Portfolios ###################################
chi_firm_dec10 <- load_and_prefix("my_chi_firm_DEC10_wide.rds", "CHI_FIRM_D")
chi_firm_exp_dec10 <- load_and_prefix("my_chi_firm_exp_DEC10_wide.rds", "CHI_FIRM_EXP_D")
chi_sic_dec10 <- load_and_prefix("my_chi_sic_DEC10_wide.rds", "CHI_SIC_D")
chi_sic_exp_dec10 <- load_and_prefix("my_chi_sic_exp_DEC10_wide.rds", "CHI_SIC_EXP_D")
chi_sic2_dec10 <- load_and_prefix("my_chi_sic2_DEC10_wide.rds", "CHI_SIC2_D")
chi_sic2_exp_dec10 <- load_and_prefix("my_chi_sic2_exp_DEC10_wide.rds", "CHI_SIC2_EXP_D")

### 5.0.5B Load Long-Short Chi Portfolios #####################################
ls_size_chi <- load_and_prefix("my_Size_CHI_LS_wide.rds", "LS_SIZE_CHI_")
ls_op_chi_sic <- load_and_prefix("my_OP_CHI_SIC_LS_wide.rds", "LS_OP_CHI_SIC_")
ls_bm_chi_sic <- load_and_prefix("my_BM_CHI_SIC_LS_wide.rds", "LS_BM_CHI_SIC_")
ls_inv_chi_sic <- load_and_prefix("my_INV_CHI_SIC_LS_wide.rds", "LS_INV_CHI_SIC_")
ls_size_within_chi <- load_and_prefix("my_Size_LS_within_CHI_wide.rds", "LS_SIZE_WITHIN_CHI_")
ls_char_op_chi_sic <- load_and_prefix("my_OP_CHI_SIC_char_LS_wide.rds", "LS_CHAR_OP_CHI_SIC_")
ls_char_bm_chi_sic <- load_and_prefix("my_BM_CHI_SIC_char_LS_wide.rds", "LS_CHAR_BM_CHI_SIC_")
ls_char_inv_chi_sic <- load_and_prefix("my_INV_CHI_SIC_char_LS_wide.rds", "LS_CHAR_INV_CHI_SIC_")

### 5.0.6 Merge Factors #######################################################
factors_all <-
  merge(
    ff5,
    SMQ,
    by = "mdate",
    all = TRUE
  )

### 5.0.7 Merge Test Assets ###################################################
test_assets <- list(

  ## 1D
  size_1d, bm1d_1d, op1d_1d, inv1d_1d, chi1d_1d,

  ## 5x5
  bm_5x5, op_5x5, inv_5x5,
  chi_firm_5x5, chi_firm_exp_5x5,
  chi_sic_5x5, chi_sic_exp_5x5,
  chi_sic2_5x5, chi_sic2_exp_5x5,
  bm_chi_sic_5x5, op_chi_sic_5x5, inv_chi_sic_5x5,

  ## 2x4x4
  bm_op_2x4x4, bm_inv_2x4x4, op_inv_2x4x4,
  op_chi_firm_2x4x4, op_chi_sic_2x4x4, op_chi_sic2_2x4x4,
  bm_chi_firm_2x4x4, bm_chi_sic_2x4x4, bm_chi_sic2_2x4x4,
  inv_chi_firm_2x4x4, inv_chi_sic_2x4x4, inv_chi_sic2_2x4x4,

  ## 1x10
  chi_firm_dec10, chi_firm_exp_dec10,
  chi_sic_dec10, chi_sic_exp_dec10,
  chi_sic2_dec10, chi_sic2_exp_dec10,

  ## Long-short chi portfolios
  ls_size_chi,
  ls_op_chi_sic,
  ls_bm_chi_sic,
  ls_inv_chi_sic,

  ## New long-short portfolios from section 4.7
  ls_size_within_chi,
  ls_char_op_chi_sic,
  ls_char_bm_chi_sic,
  ls_char_inv_chi_sic
)

test_assets_all <-
  Reduce(
    function(x, y) merge(x, y, by = "mdate", all = TRUE),
    test_assets
  )

### 5.0.8 Final Master Panel ##################################################
master_panel <-
  merge(
    factors_all,
    test_assets_all,
    by = "mdate",
    all = TRUE
  )

setkey(master_panel, mdate)
master_panel <- rf_monthly[master_panel]

factor_test_panel <-
  master_panel[
    mdate >= START_DATE &
      mdate <= END_DATE
  ]

### 5.0.9 Save Master Panel ###################################################
saveRDS(
  factor_test_panel,
  file = "master_factors_testassets.rds"
)

### 5.0A Preliminary Inspection ###############################################
cat("\n============================================================\n")
cat("MASTER FACTOR-TEST ASSET PANEL\n")
cat("============================================================\n")
print(dim(factor_test_panel))
print(range(factor_test_panel$mdate, na.rm = TRUE))
print(summary(factor_test_panel[, .(mdate, RF)]))

### 5.0B Diagnostic Checks ####################################################

cat("\n============================================================\n")
cat("DIAGNOSTIC CHECKS\n")
cat("============================================================\n")

## Check uniqueness of mdate
cat("\nUnique mdate check:\n")
print(uniqueN(factor_test_panel$mdate))
print(nrow(factor_test_panel))

## Check RF merge
cat("\nRF coverage:\n")
print(factor_test_panel[, .(
  n_months = .N,
  rf_missing = sum(is.na(RF)),
  rf_nonmissing = sum(!is.na(RF))
)])

## Check key factor coverage
cat("\nFactor coverage:\n")
factor_cols_check <- intersect(c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ"), names(factor_test_panel))
print(
  factor_test_panel[
    ,
    lapply(.SD, function(x) sum(!is.na(x))),
    .SDcols = factor_cols_check
  ]
)

## Check new 1D columns
cat("\n1D portfolio columns:\n")
print(grep("^SIZE_S[1-5]$|^BM1D_B[1-5]$|^OP1D_O[1-5]$|^INV1D_I[1-5]$|^CHI1D_C[1-5]$",
  names(factor_test_panel),
  value = TRUE
))

## Check new 5x5 chi interaction columns
cat("\n5x5 chi-interaction columns:\n")
print(grep("^BM_CHI_SIC_B[1-5]C[1-5]$|^OP_CHI_SIC_O[1-5]C[1-5]$|^INV_CHI_SIC_I[1-5]C[1-5]$",
  names(factor_test_panel),
  value = TRUE
))

## Quick missingness snapshot for new blocks
new_cols_check <- c(
  grep("^SIZE_S[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^BM1D_B[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^OP1D_O[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^INV1D_I[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^CHI1D_C[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^BM_CHI_SIC_B[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^OP_CHI_SIC_O[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  grep("^INV_CHI_SIC_I[1-5]C[1-5]$", names(factor_test_panel), value = TRUE)
)

cat("\nMissingness snapshot for new blocks:\n")
print(
  data.table(
    variable = new_cols_check,
    n_nonmissing = sapply(new_cols_check, function(v) sum(!is.na(factor_test_panel[[v]]))),
    share_nonmissing = sapply(new_cols_check, function(v) mean(!is.na(factor_test_panel[[v]])))
  )
)

cat("\nSaved: master_factors_testassets.rds\n")



master <- readRDS("master_factors_testassets.rds")
names(master)

grep("^SIZE_", names(master), value = TRUE)
grep("^BM1D_", names(master), value = TRUE)
grep("^OP1D_", names(master), value = TRUE)
grep("^INV1D_", names(master), value = TRUE)
grep("^CHI1D_", names(master), value = TRUE)

grep("^BM_CHI_SIC_", names(master), value = TRUE)
grep("^OP_CHI_SIC_", names(master), value = TRUE)
grep("^INV_CHI_SIC_", names(master), value = TRUE)
