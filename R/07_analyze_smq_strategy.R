## 7. Slow-Minus-Quick Strategy: Baseline SIC Value-Weighted ###################
## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026

### 7.0 Setup #################################################################
cat("\014")
rm(list = ls())
graphics.off()

### Libraries #################################################################
library(data.table)
library(lubridate)
library(lmtest)
library(sandwich)

### Project Paths #############################################################
source(file.path("R", "00_config.R"))

### Settings ##################################################################
START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")
NW_LAGS <- 12

### Models ####################################################################
models <- list(
  CAPM = c("MKT"),
  FF3  = c("MKT", "SMB", "HML"),
  FF5  = c("MKT", "SMB", "HML", "RMW", "CMA")
)

### 7.1 Load Data #############################################################
factor_test_panel <- readRDS("master_factors_testassets.rds")
setDT(factor_test_panel)
factor_test_panel[, mdate := as.Date(mdate)]

data <- readRDS("crsp_compustat_with_adj_cost.rds")
setDT(data)
data[, mdate := as.Date(mdate)]

### Restrict Factor Panel Window ##############################################
factor_test_panel <- factor_test_panel[
  mdate >= START_DATE &
    mdate <= END_DATE
]

### Calendar Fields ###########################################################
data[, year := year(mdate)]
data[, month := month(mdate)]
data[, ffyear := fifelse(month >= 7, year, year - 1)]

### Helpers ###################################################################
star_p <- function(p) {
  fifelse(
    p < 0.01, "***",
    fifelse(
      p < 0.05, "**",
      fifelse(p < 0.10, "*", "")
    )
  )
}

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
    pval     = as.numeric(ct[, 4])
  )
  out[, stars := star_p(pval)]
  out[]
}

mean_nw <- function(y, lags = NW_LAGS) {
  dt <- data.table(y = y)
  dt <- dt[complete.cases(dt)]

  if (nrow(dt) < 60) {
    return(NULL)
  }

  fit <- lm(y ~ 1, data = dt)
  nw <- NeweyWest(fit, lag = lags, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov. = nw)

  est <- as.numeric(ct[1, 1])
  tst <- as.numeric(ct[1, 3])
  pvl <- as.numeric(ct[1, 4])

  data.table(
    estimate = est,
    tstat    = tst,
    pval     = pvl,
    stars    = star_p(pvl)
  )
}

format_cell <- function(est, tstat, stars, digits_est = 2, digits_t = 2, scale = 100) {
  est_s <- scale * est
  sprintf(
    paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
    est_s, stars, tstat
  )
}

parse_dec10_num <- function(col, prefix) {
  x <- sub(paste0("^", prefix), "", col)
  x <- gsub("[^0-9]", "", x)
  if (x == "") {
    return(NA_integer_)
  }
  as.integer(x)
}

get_sorted_decile_cols <- function(panel, prefix) {
  cols <- grep(paste0("^", prefix), names(panel), value = TRUE)
  if (length(cols) == 0) stop("No decile columns found for prefix: ", prefix)

  nums <- vapply(cols, parse_dec10_num, integer(1), prefix = prefix)
  ok <- !is.na(nums) & nums >= 1 & nums <= 10

  cols <- cols[ok]
  nums <- nums[ok]

  ord <- order(nums)
  cols[ord]
}

cum_safe <- function(x) {
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- cumprod(1 + x[ok]) - 1
  out
}

nw_tstat_mean <- function(x, lags = NW_LAGS) {
  x <- x[!is.na(x)]
  if (length(x) < 24) {
    return(NA_real_)
  }

  fit <- lm(x ~ 1)
  nw <- NeweyWest(fit, lag = lags, prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov. = nw)

  as.numeric(ct[1, "t value"])
}

ff_30_70_breaks <- function(x) {
  qs <- quantile(x, probs = c(0.3, 0.7), na.rm = TRUE, type = 1)
  c(low = as.numeric(qs[1]), high = as.numeric(qs[2]))
}

### 7.2 Decile Table with Baseline SIC SMQ ####################################
build_decile_table_with_smq <- function(panel, dec_prefix, smq_col, lags = NW_LAGS,
                                        digits_est = 2, digits_t = 2) {
  stopifnot("RF" %in% names(panel))
  stopifnot(smq_col %in% names(panel))

  dec_cols <- get_sorted_decile_cols(panel, dec_prefix)
  col_labels <- c(paste0("P", 1:10), "SMQ")

  get_excess <- function(col) {
    panel[[col]] - panel[["RF"]]
  }

  mean_cells <- character(length(col_labels))
  vol_cells <- character(length(col_labels))
  sharpe_cells <- character(length(col_labels))

  for (j in seq_along(dec_cols)) {
    y <- get_excess(dec_cols[j])
    m <- mean_nw(y, lags)
    sd_m <- sd(y, na.rm = TRUE)

    mean_cells[j] <- if (is.null(m)) {
      ""
    } else {
      format_cell(m$estimate, m$tstat, m$stars, digits_est, digits_t, scale = 100)
    }

    vol_cells[j] <- sprintf("%.2f", 100 * sd_m)
    sharpe_cells[j] <- sprintf("%.2f", (mean(y, na.rm = TRUE) * 12) / (sd_m * sqrt(12)))
  }

  y_smq <- panel[[smq_col]]
  m_smq <- mean_nw(y_smq, lags)
  sd_smq <- sd(y_smq, na.rm = TRUE)

  mean_cells[11] <- if (is.null(m_smq)) {
    ""
  } else {
    format_cell(m_smq$estimate, m_smq$tstat, m_smq$stars, digits_est, digits_t, scale = 100)
  }
  vol_cells[11] <- sprintf("%.2f", 100 * sd_smq)
  sharpe_cells[11] <- sprintf("%.2f", (mean(y_smq, na.rm = TRUE) * 12) / (sd_smq * sqrt(12)))

  alpha_row <- function(rhs_vars) {
    out <- character(length(col_labels))
    X <- panel[, ..rhs_vars]

    for (j in seq_along(dec_cols)) {
      y <- get_excess(dec_cols[j])
      reg <- run_ts_reg_nw(y, X, lags)

      out[j] <- if (is.null(reg)) {
        ""
      } else {
        a <- reg[term == "(Intercept)"][1]
        format_cell(a$estimate, a$tstat, a$stars, digits_est, digits_t, scale = 100)
      }
    }

    reg_smq <- run_ts_reg_nw(y_smq, X, lags)
    out[11] <- if (is.null(reg_smq)) {
      ""
    } else {
      a <- reg_smq[term == "(Intercept)"][1]
      format_cell(a$estimate, a$tstat, a$stars, digits_est, digits_t, scale = 100)
    }

    out
  }

  capm_cells <- alpha_row(models$CAPM)
  ff3_cells <- alpha_row(models$FF3)
  ff5_cells <- alpha_row(models$FF5)

  tab <- rbind(
    Excess_return = mean_cells,
    Volatility    = vol_cells,
    Sharpe_ratio  = sharpe_cells,
    CAPM_alpha    = capm_cells,
    FF3_alpha     = ff3_cells,
    FF5_alpha     = ff5_cells
  )

  colnames(tab) <- col_labels
  tab
}

decile_tables <- list()

decile_tables[["CHI_SIC_DEC10_BASELINE"]] <-
  build_decile_table_with_smq(
    panel      = factor_test_panel,
    dec_prefix = "CHI_SIC_DD",
    smq_col    = "SMQ_sic_VW",
    lags       = NW_LAGS
  )

saveRDS(decile_tables, "table_bab_style_dec10_chi_sic_baseline.rds")

for (nm in names(decile_tables)) {
  cat("\n\n############################################################\n")
  cat("BASELINE SIC DECILE PORTFOLIOS AND SMQ:", nm, "\n")
  cat("############################################################\n\n")
  print(decile_tables[[nm]])
}

### 7.2B FF5 Alpha and Factor Loadings Table ##################################
build_ff5_loading_table <- function(panel, dec_prefix, smq_col, lags = NW_LAGS,
                                    digits_est = 2, digits_t = 2) {
  stopifnot("RF" %in% names(panel))
  stopifnot(smq_col %in% names(panel))

  rhs_vars <- models$FF5
  stopifnot(all(rhs_vars %in% names(panel)))

  dec_cols <- get_sorted_decile_cols(panel, dec_prefix)
  col_labels <- c(paste0("P", 1:10), "SMQ")

  get_excess <- function(col) {
    panel[[col]] - panel[["RF"]]
  }

  extract_row <- function(term_name) {
    out <- character(length(col_labels))
    X <- panel[, ..rhs_vars]

    for (j in seq_along(dec_cols)) {
      y <- get_excess(dec_cols[j])
      reg <- run_ts_reg_nw(y, X, lags)

      out[j] <- if (is.null(reg)) {
        ""
      } else {
        x <- reg[term == term_name][1]
        format_cell(x$estimate, x$tstat, x$stars,
          digits_est = digits_est,
          digits_t   = digits_t,
          scale      = ifelse(term_name == "(Intercept)", 100, 1)
        )
      }
    }

    reg_smq <- run_ts_reg_nw(panel[[smq_col]], X, lags)
    out[11] <- if (is.null(reg_smq)) {
      ""
    } else {
      x <- reg_smq[term == term_name][1]
      format_cell(x$estimate, x$tstat, x$stars,
        digits_est = digits_est,
        digits_t   = digits_t,
        scale      = ifelse(term_name == "(Intercept)", 100, 1)
      )
    }

    out
  }

  tab <- rbind(
    FF5_alpha = extract_row("(Intercept)"),
    MKT       = extract_row("MKT"),
    SMB       = extract_row("SMB"),
    HML       = extract_row("HML"),
    RMW       = extract_row("RMW"),
    CMA       = extract_row("CMA")
  )

  colnames(tab) <- col_labels
  tab
}

ff5_loading_tables <- list()

ff5_loading_tables[["CHI_SIC_DEC10_BASELINE"]] <-
  build_ff5_loading_table(
    panel      = factor_test_panel,
    dec_prefix = "CHI_SIC_DD",
    smq_col    = "SMQ_sic_VW",
    lags       = NW_LAGS
  )

saveRDS(ff5_loading_tables, "table_ff5_alpha_loadings_chi_sic_baseline.rds")

for (nm in names(ff5_loading_tables)) {
  cat("\n\n############################################################\n")
  cat("FF5 ALPHA AND FACTOR LOADINGS:", nm, "\n")
  cat("############################################################\n\n")
  print(ff5_loading_tables[[nm]])
}

### 7.2C SIC Decile Portfolio Cell Counts #####################################
labelled_data <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(labelled_data)
labelled_data[, mdate := as.Date(mdate)]

labelled_data <- labelled_data[
  mdate >= START_DATE &
    mdate <= END_DATE
]

count_tables_deciles <- list()

tmp <- labelled_data[
  !is.na(chi_sic_dec10),
  .(N = uniqueN(permno)),
  by = .(mdate, decile = chi_sic_dec10)
][
  ,
  .(Avg_N = mean(N, na.rm = TRUE)),
  by = decile
][order(decile)]

count_tables_deciles[["chi_sic_dec10"]] <- tmp

cat("\n============================================================\n")
cat("BASELINE SIC DECILE PORTFOLIO CELL COUNTS\n")
cat("============================================================\n")
print(count_tables_deciles[["chi_sic_dec10"]])

saveRDS(count_tables_deciles, "table_decile_cell_counts_sic_baseline.rds")

### 7.3 Time-Series Plots: Baseline SIC VW ####################################
par(mfrow = c(1, 1), mar = c(3, 4, 2, 1))
plot(
  factor_test_panel$mdate,
  factor_test_panel$SMQ_sic_VW,
  type = "l",
  lwd  = 1.5,
  main = "SMQ - SIC Level (VW, Baseline)",
  ylab = "Return",
  xlab = "Time"
)

par(mfrow = c(1, 1), mar = c(3, 4, 2, 1))
plot(
  factor_test_panel$mdate,
  cum_safe(factor_test_panel$SMQ_sic_VW),
  type = "l",
  lwd  = 1.5,
  main = "Cumulative SMQ - SIC Level (VW, Baseline)",
  ylab = "Cumulative Return",
  xlab = "Time"
)

### 7.4 Rebuild June Formation Sample #########################################
june_form <- data[
  month == 6 &
    has_dec == TRUE &
    has_june == TRUE &
    mdate >= START_DATE &
    mdate <= END_DATE,
  .(
    permno,
    gvkey,
    ffyear = ffyear + 1L,
    exchcd,
    bm = bm_ff,
    op = op_clean,
    inv = inv_annual,
    ik = ik_annual,
    q = tobins_q,
    chi_sic = chi_sic
  )
]

### 7.5 NYSE Breakpoints: SIC Baseline ########################################
signals_nyse <- june_form[exchcd == 1]

chi_sic_bp <- signals_nyse[
  !is.na(chi_sic),
  {
    b <- ff_30_70_breaks(chi_sic)
    .(chi_sic_30 = b["low"], chi_sic_70 = b["high"])
  },
  by = ffyear
]

june_form <- merge(
  june_form,
  chi_sic_bp,
  by = "ffyear",
  all.x = TRUE
)

### 7.6 Assign Low / High: SIC Baseline #######################################
june_form[, chi_sic_grp := fcase(
  !is.na(chi_sic) & chi_sic <= chi_sic_30, "Low",
  !is.na(chi_sic) & chi_sic > chi_sic_70, "High",
  default = NA_character_
)]

### 7.7 Factor-Leg Counts: SIC Baseline #######################################
factor_leg_counts <- june_form[
  !is.na(chi_sic_grp),
  .(N = .N),
  by = .(ffyear, Leg = chi_sic_grp)
][
  ,
  .(
    Factor   = "SMQ_sic_baseline",
    Leg,
    Avg_N    = mean(N),
    Median_N = median(N)
  )
]

cat("\n============================================================\n")
cat("BASELINE SIC FACTOR-LEG COUNTS\n")
cat("============================================================\n")
print(factor_leg_counts)

saveRDS(factor_leg_counts, "table_factor_leg_counts_sic_baseline.rds")

### 7.8 Characteristic Spread by Factor Leg: SIC Baseline #####################
make_leg_spread_table <- function(dt, grp_var, factor_name, chi_var) {
  tmp <- dt[
    !is.na(get(grp_var)),
    .(
      chi_mean = mean(get(chi_var), na.rm = TRUE),
      bm_mean  = mean(bm, na.rm = TRUE),
      op_mean  = mean(op, na.rm = TRUE),
      inv_mean = mean(inv, na.rm = TRUE),
      ik_mean  = mean(ik, na.rm = TRUE),
      q_mean   = mean(q, na.rm = TRUE)
    ),
    by = get(grp_var)
  ]

  setnames(tmp, "get", "Leg")
  tmp[, Factor := factor_name]

  spread <- tmp[
    ,
    .(
      Factor   = factor_name,
      Leg      = "High-Low",
      chi_mean = chi_mean[Leg == "High"] - chi_mean[Leg == "Low"],
      bm_mean  = bm_mean[Leg == "High"] - bm_mean[Leg == "Low"],
      op_mean  = op_mean[Leg == "High"] - op_mean[Leg == "Low"],
      inv_mean = inv_mean[Leg == "High"] - inv_mean[Leg == "Low"],
      ik_mean  = ik_mean[Leg == "High"] - ik_mean[Leg == "Low"],
      q_mean   = q_mean[Leg == "High"] - q_mean[Leg == "Low"]
    )
  ]

  rbind(
    tmp[, .(Factor, Leg, chi_mean, bm_mean, op_mean, inv_mean, ik_mean, q_mean)],
    spread
  )
}

leg_spread_table <- make_leg_spread_table(
  june_form,
  "chi_sic_grp",
  "SMQ_sic_baseline",
  "chi_sic"
)

cat("\n============================================================\n")
cat("BASELINE SIC CHARACTERISTIC SPREADS BY FACTOR LEG\n")
cat("============================================================\n")
print(leg_spread_table)

saveRDS(
  leg_spread_table,
  "table_factor_leg_characteristic_spreads_sic_baseline.rds"
)

### 7.9 Merge Leg Labels into Monthly Panel ###################################
leg_labels <- june_form[
  ,
  .(permno, ffyear, chi_sic_grp)
]

monthly <- merge(
  data,
  leg_labels,
  by = c("permno", "ffyear"),
  all.x = TRUE
)

setorder(monthly, permno, mdate)

monthly[, me_ff_w := shift(me_clean), by = permno]
monthly[me_ff_w <= 0, me_ff_w := NA_real_]

### 7.10 Factor-Leg Excess Return Statistics: SIC Baseline ####################
make_leg_return_stats <- function(dt, grp_var, factor_name) {
  leg_rets <- dt[
    mdate >= START_DATE &
      mdate <= END_DATE &
      !is.na(get(grp_var)) &
      !is.na(me_ff_w) &
      me_ff_w > 0 &
      !is.na(retadj) &
      !is.na(RF),
    .(
      vwret = weighted.mean(retadj, me_ff_w, na.rm = TRUE),
      RF    = mean(RF, na.rm = TRUE)
    ),
    by = .(mdate, Leg = get(grp_var))
  ]

  leg_rets[, exret := vwret - RF]

  spread <- leg_rets[
    ,
    .(exret = exret[Leg == "High"] - exret[Leg == "Low"]),
    by = mdate
  ]
  spread[, Leg := "High-Low"]

  all_rets <- rbind(
    leg_rets[, .(mdate, Leg, exret)],
    spread[, .(mdate, Leg, exret)],
    fill = TRUE
  )

  out <- all_rets[
    ,
    .(
      Mean_excess = mean(exret, na.rm = TRUE),
      Tstat = nw_tstat_mean(exret, NW_LAGS),
      Vol = sd(exret, na.rm = TRUE) * sqrt(12),
      Sharpe = ifelse(
        sd(exret, na.rm = TRUE) > 0,
        mean(exret, na.rm = TRUE) / sd(exret, na.rm = TRUE) * sqrt(12),
        NA_real_
      )
    ),
    by = Leg
  ]

  out[, Factor := factor_name]
  out[, Mean_excess := 100 * Mean_excess]
  out[, Vol := 100 * Vol]

  out[, .(Factor, Leg, Mean_excess, Tstat, Vol, Sharpe)]
}

factor_leg_return_stats <- make_leg_return_stats(
  monthly,
  "chi_sic_grp",
  "SMQ_sic_baseline"
)

cat("\n============================================================\n")
cat("BASELINE SIC FACTOR-LEG EXCESS RETURN STATISTICS\n")
cat("============================================================\n")
print(factor_leg_return_stats)

saveRDS(
  factor_leg_return_stats,
  "table_factor_leg_excess_return_stats_sic_baseline.rds"
)

### 7.11 High-Low Spread Series: SIC Baseline #################################
make_leg_spread_series <- function(dt, grp_var, factor_name) {
  leg_rets <- dt[
    mdate >= START_DATE &
      mdate <= END_DATE &
      !is.na(get(grp_var)) &
      !is.na(me_ff_w) &
      me_ff_w > 0 &
      !is.na(retadj) &
      !is.na(RF),
    .(
      vwret = weighted.mean(retadj, me_ff_w, na.rm = TRUE),
      RF    = mean(RF, na.rm = TRUE)
    ),
    by = .(mdate, Leg = get(grp_var))
  ]

  leg_rets[, exret := vwret - RF]

  spread <- leg_rets[
    ,
    .(exret = exret[Leg == "High"] - exret[Leg == "Low"]),
    by = mdate
  ]

  spread[, Factor := factor_name]
  spread
}

spread_sic <- make_leg_spread_series(
  monthly,
  "chi_sic_grp",
  "SMQ_sic_baseline"
)

saveRDS(
  spread_sic,
  "spread_series_smq_sic_high_low_baseline.rds"
)

### 7.12 Distribution Plots: SIC Baseline #####################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
hist(
  100 * spread_sic$exret,
  breaks = 30,
  main   = "High-Low Excess Return\nSMQ_sic Baseline",
  xlab   = "Monthly excess return (%)"
)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(
  density(100 * spread_sic$exret, na.rm = TRUE),
  main = "Kernel Density\nSMQ_sic Baseline",
  xlab = "Monthly excess return (%)"
)
abline(v = mean(100 * spread_sic$exret, na.rm = TRUE), lty = 2)

### 7.12B Joint Plot: Cumulative Return and Kernel Density ####################
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

## Left panel: cumulative SMQ return
plot(
  factor_test_panel$mdate,
  cum_safe(factor_test_panel$SMQ_sic_VW),
  type = "l",
  lwd  = 1.5,
  main = "Cumulative SMQ - SIC Level (VW, Baseline)",
  ylab = "Cumulative Return",
  xlab = "Time"
)
abline(h = 0, lty = 2)

## Right panel: kernel density of monthly SMQ returns
plot(
  density(100 * factor_test_panel$SMQ_sic_VW, na.rm = TRUE),
  main = "Kernel Density of SMQ - SIC Level",
  xlab = "Monthly return (%)",
  lwd = 1.5
)
abline(v = mean(100 * factor_test_panel$SMQ_sic_VW, na.rm = TRUE), lty = 2)

par(mfrow = c(1, 1))

### 7.12C Normality Diagnostics: SIC Baseline #################################
make_dist_summary <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 24) {
    return(NULL)
  }

  mu <- mean(x)
  sd_x <- sd(x)

  skew <- mean((x - mu)^3) / (sd_x^3)
  kurt <- mean((x - mu)^4) / (sd_x^4)
  ex_kurt <- kurt - 3

  jb <- (length(x) / 6) * skew^2 + (length(x) / 24) * (ex_kurt)^2
  jb_p <- 1 - pchisq(jb, df = 2)

  out <- data.table(
    N               = length(x),
    Mean            = mu,
    SD              = sd_x,
    Skewness        = skew,
    Kurtosis        = kurt,
    Excess_kurtosis = ex_kurt,
    Min             = min(x),
    P01             = as.numeric(quantile(x, 0.01)),
    P05             = as.numeric(quantile(x, 0.05)),
    P25             = as.numeric(quantile(x, 0.25)),
    Median          = as.numeric(quantile(x, 0.50)),
    P75             = as.numeric(quantile(x, 0.75)),
    P95             = as.numeric(quantile(x, 0.95)),
    P99             = as.numeric(quantile(x, 0.99)),
    Max             = max(x),
    JB_stat         = jb,
    JB_pval         = jb_p
  )

  out
}


smq_raw <- factor_test_panel$SMQ_sic_VW
smq_raw <- smq_raw[!is.na(smq_raw)]

smq_z <- (smq_raw - mean(smq_raw)) / sd(smq_raw)

smq_dist_summary_raw <- make_dist_summary(smq_raw)
smq_dist_summary_z <- make_dist_summary(smq_z)

cat("\n============================================================\n")
cat("SMQ SIC BASELINE DISTRIBUTION SUMMARY: RAW RETURNS\n")
cat("============================================================\n")
print(round(smq_dist_summary_raw, 4))

cat("\n============================================================\n")
cat("SMQ SIC BASELINE DISTRIBUTION SUMMARY: STANDARDIZED RETURNS\n")
cat("============================================================\n")
print(round(smq_dist_summary_z, 4))

saveRDS(
  list(
    raw_returns          = smq_dist_summary_raw,
    standardized_returns = smq_dist_summary_z
  ),
  "table_smq_sic_baseline_normality_diagnostics.rds"
)

### 7.12D Density vs Normal: Raw Returns ######################################
smq_mean_pct <- mean(100 * smq_raw, na.rm = TRUE)
smq_sd_pct <- sd(100 * smq_raw, na.rm = TRUE)

dens_raw <- density(100 * smq_raw, na.rm = TRUE)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(
  dens_raw,
  lwd  = 1.5,
  main = "SMQ_sic Baseline: Kernel Density vs Fitted Normal",
  xlab = "Monthly return (%)"
)
curve(
  dnorm(x, mean = smq_mean_pct, sd = smq_sd_pct),
  add = TRUE,
  lty = 2,
  lwd = 1.5
)
abline(v = smq_mean_pct, lty = 3)
legend(
  "topright",
  legend = c("Kernel density", "Fitted normal", "Sample mean"),
  lty    = c(1, 2, 3),
  lwd    = c(1.5, 1.5, 1),
  bty    = "n"
)

### 7.12E Density vs Standard Normal: Standardized Returns ####################
dens_z <- density(smq_z, na.rm = TRUE)

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(
  dens_z,
  lwd  = 1.5,
  main = "Standardized SMQ_sic: Kernel Density vs Standard Normal",
  xlab = "Standardized monthly return"
)
curve(
  dnorm(x, mean = 0, sd = 1),
  add = TRUE,
  lty = 2,
  lwd = 1.5
)
abline(v = 0, lty = 3)
legend(
  "topright",
  legend = c("Kernel density", "Standard normal", "Zero"),
  lty    = c(1, 2, 3),
  lwd    = c(1.5, 1.5, 1),
  bty    = "n"
)

### 7.12F QQ Plot: Standardized Returns #######################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
qqnorm(
  smq_z,
  main = "QQ Plot: Standardized SMQ_sic Baseline",
  ylab = "Sample quantiles"
)
qqline(smq_z, lty = 2)

### 7.12G Histogram vs Fitted Normal ##########################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
hist(
  100 * smq_raw,
  breaks = 30,
  freq   = FALSE,
  main   = "SMQ_sic Baseline: Histogram vs Fitted Normal",
  xlab   = "Monthly return (%)"
)
curve(
  dnorm(x, mean = smq_mean_pct, sd = smq_sd_pct),
  add = TRUE,
  lty = 2,
  lwd = 1.5
)
abline(v = smq_mean_pct, lty = 3)

### 7.12H Joint Plot: Raw Density and QQ Plot #################################
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(
  dens_raw,
  lwd  = 1.5,
  main = "Kernel Density vs Fitted Normal",
  xlab = "Monthly return (%)"
)
curve(
  dnorm(x, mean = smq_mean_pct, sd = smq_sd_pct),
  add = TRUE,
  lty = 2,
  lwd = 1.5
)
abline(v = smq_mean_pct, lty = 3)

qqnorm(
  smq_z,
  main = "QQ Plot: Standardized SMQ",
  ylab = "Sample quantiles"
)
qqline(smq_z, lty = 2)

par(mfrow = c(1, 1))

### 7.13 Average Chi Level by Decile ##########################################
labelled_data <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(labelled_data)
labelled_data[, mdate := as.Date(mdate)]

labelled_data <- labelled_data[
  mdate >= START_DATE &
    mdate <= END_DATE
]

chi_by_decile <- labelled_data[
  !is.na(chi_sic_dec10) &
    !is.na(chi_sic),
  .(
    Avg_chi    = mean(chi_sic, na.rm = TRUE),
    Median_chi = median(chi_sic, na.rm = TRUE),
    SD_chi     = sd(chi_sic, na.rm = TRUE),
    N_obs      = .N
  ),
  by = chi_sic_dec10
][order(chi_sic_dec10)]

cat("\n============================================================\n")
cat("AVERAGE CHI LEVEL BY SIC DECILE\n")
cat("============================================================\n")
print(chi_by_decile)

saveRDS(chi_by_decile, "table_avg_chi_by_sic_decile.rds")

### 7.14 Industry Concentration (HHI) #########################################
## Uses SIC2 as the baseline industry grouping for concentration diagnostics.
## HHI = sum of squared portfolio industry shares within each decile-month,
## then averaged across time.

if (!"sic2" %in% names(labelled_data)) {
  labelled_data[, sic2 := fifelse(!is.na(sic), sic %/% 100L, NA_integer_)]
}

decile_industry_weights <- labelled_data[
  !is.na(chi_sic_dec10) &
    !is.na(sic2),
  .(N_firms = uniqueN(permno)),
  by = .(mdate, chi_sic_dec10, sic2)
]

decile_industry_weights[
  ,
  total_firms := sum(N_firms),
  by = .(mdate, chi_sic_dec10)
]

decile_industry_weights[
  ,
  share := fifelse(total_firms > 0, N_firms / total_firms, NA_real_)
]

hhi_by_month <- decile_industry_weights[
  ,
  .(HHI = sum(share^2, na.rm = TRUE)),
  by = .(mdate, chi_sic_dec10)
]

hhi_by_decile <- hhi_by_month[
  ,
  .(
    Avg_HHI    = mean(HHI, na.rm = TRUE),
    Median_HHI = median(HHI, na.rm = TRUE),
    SD_HHI     = sd(HHI, na.rm = TRUE),
    N_months   = .N
  ),
  by = chi_sic_dec10
][order(chi_sic_dec10)]

cat("\n============================================================\n")
cat("INDUSTRY CONCENTRATION (HHI) BY SIC DECILE\n")
cat("============================================================\n")
print(hhi_by_decile)

saveRDS(hhi_by_decile, "table_industry_hhi_by_sic_decile.rds")

### 7.15 Dominant Industries Within Each Decile ################################
## Show which SIC2 industries dominate each decile on average across months.

industry_shares_avg <- decile_industry_weights[
  ,
  .(Avg_share = mean(share, na.rm = TRUE)),
  by = .(chi_sic_dec10, sic2)
][order(chi_sic_dec10, -Avg_share)]

top_industries_by_decile <- industry_shares_avg[
  ,
  head(.SD, 10),
  by = chi_sic_dec10
]

cat("\n============================================================\n")
cat("TOP 10 SIC2 INDUSTRIES BY AVERAGE SHARE WITHIN EACH DECILE\n")
cat("============================================================\n")
print(top_industries_by_decile)

saveRDS(top_industries_by_decile, "table_top_industries_by_sic_decile.rds")

### 7.15B Industry-Level Characteristics Within Each Decile ####################
## Average firm characteristics within each decile-SIC2 cell.
## Helps assess whether certain industry types dominate the extreme deciles.

industry_char_by_decile <- labelled_data[
  !is.na(chi_sic_dec10) &
    !is.na(sic2),
  .(
    Avg_chi = mean(chi_sic, na.rm = TRUE),
    Avg_bm  = mean(bm_ff, na.rm = TRUE),
    Avg_op  = mean(op_clean, na.rm = TRUE),
    Avg_inv = mean(inv_annual, na.rm = TRUE),
    Avg_q   = mean(tobins_q, na.rm = TRUE),
    N_obs   = .N
  ),
  by = .(chi_sic_dec10, sic2)
][order(chi_sic_dec10, -N_obs)]

saveRDS(industry_char_by_decile, "table_industry_characteristics_by_sic_decile.rds")

### 7.16 Correlation Matrix of SIC Decile Portfolio Excess Returns + SMQ ######
dec_cols <- get_sorted_decile_cols(factor_test_panel, "CHI_SIC_DD")

corr_dt <- factor_test_panel[
  ,
  c("mdate", "RF", dec_cols, "SMQ_sic_VW"),
  with = FALSE
]

corr_excess <- copy(corr_dt)

for (col in dec_cols) {
  corr_excess[, (col) := get(col) - RF]
}

## SMQ_sic_VW is already a long-short factor return, so do NOT subtract RF
corr_vars <- c(dec_cols, "SMQ_sic_VW")

corr_matrix <- cor(
  corr_excess[, ..corr_vars],
  use = "complete.obs"
)

rownames(corr_matrix) <- c(paste0("P", 1:10), "SMQ")
colnames(corr_matrix) <- c(paste0("P", 1:10), "SMQ")

cat("\n============================================================\n")
cat("CORRELATION MATRIX OF SIC DECILE PORTFOLIO EXCESS RETURNS AND SMQ\n")
cat("============================================================\n")
print(round(corr_matrix, 4))

saveRDS(
  corr_matrix,
  "table_correlation_matrix_sic_decile_excess_returns_with_smq.rds"
)


### 7.17 Time-Series Plots of SIC Decile Portfolios ############################
dec_cols <- get_sorted_decile_cols(factor_test_panel, "CHI_SIC_DD")

## 7.17A Individual decile return plots
par(mfrow = c(5, 2), mar = c(3, 4, 2, 1))
for (col in dec_cols) {
  plot(
    factor_test_panel$mdate,
    factor_test_panel[[col]],
    type = "l",
    lwd  = 1.2,
    main = col,
    xlab = "Time",
    ylab = "Return"
  )
  abline(h = 0, lty = 2)
}

## 7.17B Combined decile return plot
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
matplot(
  factor_test_panel$mdate,
  as.matrix(factor_test_panel[, ..dec_cols]),
  type = "l",
  lty  = 1,
  lwd  = 1.2,
  xlab = "Time",
  ylab = "Return",
  main = "SIC Decile Portfolio Returns"
)
legend(
  "topright",
  legend = paste0("P", 1:10),
  col    = 1:10,
  lty    = 1,
  cex    = 0.8,
  bty    = "n"
)

## 7.17C Combined cumulative return plot
cum_dec_mat <- copy(factor_test_panel[, ..dec_cols])
for (j in seq_along(dec_cols)) {
  cum_dec_mat[[dec_cols[j]]] <- cum_safe(cum_dec_mat[[dec_cols[j]]])
}

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
matplot(
  factor_test_panel$mdate,
  as.matrix(cum_dec_mat),
  type = "l",
  lty  = 1,
  lwd  = 1.2,
  xlab = "Time",
  ylab = "Cumulative Return",
  main = "Cumulative SIC Decile Portfolio Returns"
)
legend(
  "topleft",
  legend = paste0("P", 1:10),
  col    = 1:10,
  lty    = 1,
  cex    = 0.8,
  bty    = "n"
)


### 7.18 Mean-Variance Comparison: FF5 vs FF5 + SMQ ############################

### 7.18A Helper Functions #####################################################
mean_monthly <- function(x) {
  mean(x, na.rm = TRUE)
}

sd_monthly <- function(x) {
  sd(x, na.rm = TRUE)
}

sharpe_annual_from_monthly <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sig <- sd(x, na.rm = TRUE)
  if (is.na(sig) || sig <= 0) {
    return(NA_real_)
  }
  sqrt(12) * (mu / sig)
}

compute_asset_stats <- function(dt, asset_names) {
  out <- rbindlist(lapply(asset_names, function(v) {
    x <- dt[[v]]
    data.table(
      Asset = v,
      Mean_excess = 100 * mean_monthly(x),
      Volatility = 100 * sd_monthly(x),
      Sharpe_ratio = sharpe_annual_from_monthly(x)
    )
  }))

  out[]
}

compute_tangency_portfolio <- function(dt, asset_names) {
  X <- dt[, ..asset_names]
  X <- X[complete.cases(X)]

  mu <- colMeans(X)
  Sigma <- cov(X)

  inv_Sigma <- solve(Sigma)
  w_unnorm <- inv_Sigma %*% mu
  w <- as.numeric(w_unnorm / sum(w_unnorm))
  names(w) <- asset_names

  port_excess <- as.numeric(as.matrix(X) %*% w)

  list(
    weights = data.table(
      Asset = asset_names,
      Weight = w
    ),
    returns = port_excess,
    summary = data.table(
      Mean_excess = 100 * mean(port_excess, na.rm = TRUE),
      Volatility = 100 * sd(port_excess, na.rm = TRUE),
      Sharpe_ratio = ifelse(
        sd(port_excess, na.rm = TRUE) > 0,
        sqrt(12) * mean(port_excess, na.rm = TRUE) / sd(port_excess, na.rm = TRUE),
        NA_real_
      )
    ),
    mu = mu,
    Sigma = Sigma
  )
}

compute_efficient_frontier <- function(dt, asset_names, n_points = 300) {
  X <- dt[, ..asset_names]
  X <- X[complete.cases(X)]

  mu <- colMeans(X)
  Sigma <- cov(X)
  inv_Sigma <- solve(Sigma)
  one <- rep(1, length(asset_names))

  A <- as.numeric(t(one) %*% inv_Sigma %*% one)
  B <- as.numeric(t(one) %*% inv_Sigma %*% mu)
  C <- as.numeric(t(mu) %*% inv_Sigma %*% mu)
  D <- A * C - B^2

  if (D <= 0) stop("Efficient frontier failed because D <= 0.")

  mu_tan <- as.numeric(
    t((inv_Sigma %*% mu) / sum(inv_Sigma %*% mu)) %*% mu
  )

  target_mu_seq <- seq(
    from = min(mu) * 0.5,
    to = max(mu_tan * 1.25, max(mu) * 1.5),
    length.out = n_points
  )

  frontier <- rbindlist(lapply(target_mu_seq, function(target_mu) {
    lambda1 <- (C - B * target_mu) / D
    lambda2 <- (A * target_mu - B) / D

    w <- inv_Sigma %*% (lambda1 * one + lambda2 * mu)
    w <- as.numeric(w)

    port_mu <- sum(w * mu)
    port_sd <- sqrt(as.numeric(t(w) %*% Sigma %*% w))

    data.table(
      Mean_excess = 100 * port_mu,
      Volatility = 100 * port_sd
    )
  }))

  frontier[]
}

print_dt_round <- function(dt, digits = 4) {
  out <- copy(dt)
  num_cols <- names(out)[sapply(out, is.numeric)]
  for (j in num_cols) {
    out[, (j) := round(get(j), digits)]
  }
  print(out)
}

### 7.18B Prepare Data #########################################################
mv_dt <- copy(factor_test_panel)

ff5_assets <- c("MKT", "SMB", "HML", "RMW", "CMA")
ff6_assets <- c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_VW")

mv_dt <- mv_dt[
  mdate >= START_DATE &
    mdate <= END_DATE
]

mv_dt <- mv_dt[
  complete.cases(mv_dt[, c(ff5_assets, "SMQ_sic_VW"), with = FALSE])
]

### 7.18C Individual Factor Statistics ########################################
ff5_factor_stats <- compute_asset_stats(mv_dt, ff5_assets)
ff6_factor_stats <- compute_asset_stats(mv_dt, ff6_assets)

cat("\n============================================================\n")
cat("FACTOR STATISTICS: FIVE-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(ff5_factor_stats)

cat("\n============================================================\n")
cat("FACTOR STATISTICS: SIX-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(ff6_factor_stats)

saveRDS(ff5_factor_stats, "table_mv_factor_stats_ff5.rds")
saveRDS(ff6_factor_stats, "table_mv_factor_stats_ff5_plus_smq.rds")

### 7.18D Tangency Portfolio: Five-Factor ######################################
tan_ff5 <- compute_tangency_portfolio(mv_dt, ff5_assets)

cat("\n============================================================\n")
cat("TANGENCY PORTFOLIO WEIGHTS: FIVE-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(tan_ff5$weights)

cat("\n============================================================\n")
cat("TANGENCY PORTFOLIO SUMMARY: FIVE-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(tan_ff5$summary)

saveRDS(tan_ff5$weights, "table_tangency_weights_ff5.rds")
saveRDS(tan_ff5$summary, "table_tangency_summary_ff5.rds")

### 7.18E Tangency Portfolio: Six-Factor #######################################
tan_ff6 <- compute_tangency_portfolio(mv_dt, ff6_assets)

cat("\n============================================================\n")
cat("TANGENCY PORTFOLIO WEIGHTS: SIX-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(tan_ff6$weights)

cat("\n============================================================\n")
cat("TANGENCY PORTFOLIO SUMMARY: SIX-FACTOR MODEL\n")
cat("============================================================\n")
print_dt_round(tan_ff6$summary)

saveRDS(tan_ff6$weights, "table_tangency_weights_ff5_plus_smq.rds")
saveRDS(tan_ff6$summary, "table_tangency_summary_ff5_plus_smq.rds")

### 7.18F Consolidated Mean-Variance Table #####################################
ff5_factor_table <- copy(ff5_factor_stats)
ff5_factor_table[, Model := "Five-Factor"]
ff5_factor_table[, Type := "Factor"]
ff5_factor_table[, Weight := NA_real_]
setcolorder(ff5_factor_table, c("Model", "Type", "Asset", "Weight", "Mean_excess", "Volatility", "Sharpe_ratio"))

ff6_factor_table <- copy(ff6_factor_stats)
ff6_factor_table[, Model := "Six-Factor"]
ff6_factor_table[, Type := "Factor"]
ff6_factor_table[, Weight := NA_real_]
setcolorder(ff6_factor_table, c("Model", "Type", "Asset", "Weight", "Mean_excess", "Volatility", "Sharpe_ratio"))

ff5_tan_table <- copy(tan_ff5$weights)
ff5_tan_table[, Mean_excess := tan_ff5$summary$Mean_excess]
ff5_tan_table[, Volatility := tan_ff5$summary$Volatility]
ff5_tan_table[, Sharpe_ratio := tan_ff5$summary$Sharpe_ratio]
ff5_tan_table[, Model := "Five-Factor"]
ff5_tan_table[, Type := "Tangency weight"]
setcolorder(ff5_tan_table, c("Model", "Type", "Asset", "Weight", "Mean_excess", "Volatility", "Sharpe_ratio"))

ff6_tan_table <- copy(tan_ff6$weights)
ff6_tan_table[, Mean_excess := tan_ff6$summary$Mean_excess]
ff6_tan_table[, Volatility := tan_ff6$summary$Volatility]
ff6_tan_table[, Sharpe_ratio := tan_ff6$summary$Sharpe_ratio]
ff6_tan_table[, Model := "Six-Factor"]
ff6_tan_table[, Type := "Tangency weight"]
setcolorder(ff6_tan_table, c("Model", "Type", "Asset", "Weight", "Mean_excess", "Volatility", "Sharpe_ratio"))
tangency_summary_table <- rbindlist(list(
  data.table(
    Model = "Five-Factor",
    Type = "Tangency summary",
    Asset = "Tangency portfolio",
    Weight = NA_real_,
    Mean_excess = tan_ff5$summary$Mean_excess,
    Volatility = tan_ff5$summary$Volatility,
    Sharpe_ratio = tan_ff5$summary$Sharpe_ratio
  ),
  data.table(
    Model = "Six-Factor",
    Type = "Tangency summary",
    Asset = "Tangency portfolio",
    Weight = NA_real_,
    Mean_excess = tan_ff6$summary$Mean_excess,
    Volatility = tan_ff6$summary$Volatility,
    Sharpe_ratio = tan_ff6$summary$Sharpe_ratio
  )
), fill = TRUE)

mv_summary_table <- rbindlist(list(
  ff5_factor_table,
  ff6_factor_table,
  ff5_tan_table,
  ff6_tan_table,
  tangency_summary_table
), fill = TRUE)

model_order <- c("Five-Factor", "Six-Factor")
type_order <- c("Factor", "Tangency weight", "Tangency summary")

mv_summary_table[, Model := factor(Model, levels = model_order)]
mv_summary_table[, Type := factor(Type, levels = type_order)]

setorder(mv_summary_table, Model, Type, Asset)

cat("\n============================================================\n")
cat("CONSOLIDATED MEAN-VARIANCE TABLE\n")
cat("============================================================\n")
print_dt_round(mv_summary_table)

saveRDS(mv_summary_table, "table_mean_variance_summary_ff5_vs_ff5_plus_smq.rds")

### 7.18G Frontier Computation #################################################
frontier_ff5 <- compute_efficient_frontier(mv_dt, ff5_assets, n_points = 300)
frontier_ff6 <- compute_efficient_frontier(mv_dt, ff6_assets, n_points = 300)

saveRDS(frontier_ff5, "frontier_ff5.rds")
saveRDS(frontier_ff6, "frontier_ff5_plus_smq.rds")

### 7.18H Mean-Variance Plotting Functions #####################################
prepare_factor_points <- function(ff5_factor_stats, ff6_factor_stats) {
  ff5_pts <- copy(ff5_factor_stats)
  ff6_pts <- copy(ff6_factor_stats)

  ff5_pts[, Label := fifelse(
    Asset == "MKT", "MKT",
    fifelse(
      Asset == "SMB", "SMB",
      fifelse(
        Asset == "HML", "HML",
        fifelse(
          Asset == "RMW", "RMW",
          fifelse(Asset == "CMA", "CMA", Asset)
        )
      )
    )
  )]

  ff6_pts[, Label := fifelse(
    Asset == "SMQ_sic_VW", "SMQ",
    fifelse(
      Asset == "MKT", "MKT",
      fifelse(
        Asset == "SMB", "SMB",
        fifelse(
          Asset == "HML", "HML",
          fifelse(
            Asset == "RMW", "RMW",
            fifelse(Asset == "CMA", "CMA", Asset)
          )
        )
      )
    )
  )]

  list(ff5 = ff5_pts, ff6 = ff6_pts)
}

get_mv_limits <- function(frontiers, tangencies, factor_sets) {
  x_max <- max(
    unlist(lapply(frontiers, function(x) x$Volatility)),
    unlist(lapply(tangencies, function(x) x$summary$Volatility)),
    unlist(lapply(factor_sets, function(x) x$Volatility)),
    na.rm = TRUE
  ) * 1.15

  y_max <- max(
    unlist(lapply(frontiers, function(x) x$Mean_excess)),
    unlist(lapply(tangencies, function(x) x$summary$Mean_excess)),
    unlist(lapply(factor_sets, function(x) x$Mean_excess)),
    na.rm = TRUE
  ) * 1.15

  y_min <- min(
    unlist(lapply(frontiers, function(x) x$Mean_excess)),
    unlist(lapply(factor_sets, function(x) x$Mean_excess)),
    0,
    na.rm = TRUE
  )

  list(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max)
  )
}

add_factor_labels <- function(dt, pos_map = NULL, cex = 0.8, xoff = 0, yoff = 0) {
  for (j in seq_len(nrow(dt))) {
    lab <- dt$Label[j]
    pos_j <- 4

    if (!is.null(pos_map) && lab %in% names(pos_map)) {
      pos_j <- pos_map[[lab]]
    }

    text(
      x = dt$Volatility[j] + xoff,
      y = dt$Mean_excess[j] + yoff,
      labels = lab,
      pos = pos_j,
      cex = cex
    )
  }
}

plot_mv_ff5 <- function(frontier_ff5, tangency_ff5, ff5_factor_stats,
                        xlim, ylim) {
  plot(
    frontier_ff5$Volatility,
    frontier_ff5$Mean_excess,
    type = "l",
    lwd = 2.7,
    lty = 1,
    xlim = xlim,
    ylim = ylim,
    xlab = "Monthly volatility (%)",
    ylab = "Monthly mean excess return (%)",
    main = "Mean-Variance Space: Five-Factor Model"
  )

  points(
    ff5_factor_stats$Volatility,
    ff5_factor_stats$Mean_excess,
    pch = 16,
    cex = 1.15
  )

  add_factor_labels(
    ff5_factor_stats,
    pos_map = list(MKT = 4, SMB = 4, HML = 4, RMW = 4, CMA = 4),
    cex = 0.8
  )

  points(
    tangency_ff5$summary$Volatility,
    tangency_ff5$summary$Mean_excess,
    pch = 17,
    cex = 1.55
  )

  text(
    tangency_ff5$summary$Volatility - 0.12,
    tangency_ff5$summary$Mean_excess + 0.03,
    labels = "Tangency portfolio",
    pos = 4,
    cex = 0.8
  )

  legend(
    "topleft",
    legend = c("Five-Factor Frontier"),
    lty = c(1),
    lwd = c(2.7),
    bty = "n",
    cex = 0.9
  )
}

plot_mv_ff6 <- function(frontier_ff6, tangency_ff6, ff6_factor_stats,
                        xlim, ylim) {
  plot(
    frontier_ff6$Volatility,
    frontier_ff6$Mean_excess,
    type = "l",
    lwd = 2.7,
    lty = 1,
    xlim = xlim,
    ylim = ylim,
    xlab = "Monthly volatility (%)",
    ylab = "Monthly mean excess return (%)",
    main = "Mean-Variance Space: Six-Factor Model"
  )

  ff6_non_smq <- ff6_factor_stats[Label != "SMQ"]
  ff6_smq <- ff6_factor_stats[Label == "SMQ"]

  points(
    ff6_non_smq$Volatility,
    ff6_non_smq$Mean_excess,
    pch = 16,
    cex = 1.15
  )

  add_factor_labels(
    ff6_non_smq,
    cex = 0.8
  )

  if (nrow(ff6_smq) == 1) {
    points(
      ff6_smq$Volatility,
      ff6_smq$Mean_excess,
      pch = 8,
      cex = 1.6
    )

    text(
      ff6_smq$Volatility + 0.05,
      ff6_smq$Mean_excess,
      labels = "SMQ",
      pos = 4,
      cex = 0.8
    )
  }

  points(
    tangency_ff6$summary$Volatility,
    tangency_ff6$summary$Mean_excess,
    pch = 17,
    cex = 1.55
  )

  text(
    tangency_ff6$summary$Volatility + 0.06,
    tangency_ff6$summary$Mean_excess,
    labels = "Tangency portfolio",
    pos = 4,
    cex = 0.8
  )

  legend(
    "topleft",
    legend = c("Six-Factor Frontier"),
    lty = c(1),
    lwd = c(2.7),
    bty = "n",
    cex = 0.9
  )
}

plot_mv_both <- function(frontier_ff5, frontier_ff6,
                         tangency_ff5, tangency_ff6,
                         ff5_factor_stats, ff6_factor_stats,
                         xlim, ylim) {
  plot(
    frontier_ff5$Volatility,
    frontier_ff5$Mean_excess,
    type = "l",
    lwd = 2.7,
    lty = 1,
    col = "black",
    xlim = xlim,
    ylim = ylim,
    xlab = "Monthly volatility (%)",
    ylab = "Monthly mean excess return (%)",
    main = "Mean-Variance Space: Five-Factor vs Six-Factor Model"
  )

  lines(
    frontier_ff6$Volatility,
    frontier_ff6$Mean_excess,
    lwd = 2.7,
    lty = 1,
    col = "gray45"
  )

  points(
    ff5_factor_stats$Volatility,
    ff5_factor_stats$Mean_excess,
    pch = 1,
    cex = 1.2
  )

  add_factor_labels(
    ff5_factor_stats,
    pos_map = list(MKT = 4, SMB = 4, HML = 4, RMW = 4, CMA = 4),
    cex = 0.8
  )

  smq_point <- ff6_factor_stats[Label == "SMQ"]
  if (nrow(smq_point) == 1) {
    points(
      smq_point$Volatility,
      smq_point$Mean_excess,
      pch = 8,
      cex = 1.6
    )

    text(
      smq_point$Volatility + 0.05,
      smq_point$Mean_excess,
      labels = "SMQ",
      pos = 4,
      cex = 0.8
    )
  }

  points(
    tangency_ff5$summary$Volatility,
    tangency_ff5$summary$Mean_excess,
    pch = 17,
    cex = 1.55
  )

  points(
    tangency_ff6$summary$Volatility,
    tangency_ff6$summary$Mean_excess,
    pch = 17,
    cex = 1.55
  )

  text(
    tangency_ff5$summary$Volatility,
    tangency_ff5$summary$Mean_excess + 0.06,
    labels = "Five-Factor Tangency Portfolio",
    pos = 4,
    cex = 0.8
  )

  text(
    tangency_ff6$summary$Volatility + 0.08,
    tangency_ff6$summary$Mean_excess,
    labels = "Six-Factor Tangency Portfolio",
    pos = 4,
    cex = 0.8
  )

  legend(
    "topleft",
    legend = c(
      "Five-Factor Frontier",
      "Six-Factor Frontier"
    ),
    lty = c(1, 1),
    lwd = c(2.7, 2.7),
    col = c("black", "gray45"),
    bty = "n",
    cex = 0.9
  )
}

### 7.18I Prepare labels and common axis limits ################################
factor_points <- prepare_factor_points(ff5_factor_stats, ff6_factor_stats)

mv_limits <- get_mv_limits(
  frontiers = list(frontier_ff5, frontier_ff6),
  tangencies = list(tan_ff5, tan_ff6),
  factor_sets = list(factor_points$ff5, factor_points$ff6)
)

### 7.18J Plot 1: Five-Factor ##################################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot_mv_ff5(
  frontier_ff5 = frontier_ff5,
  tangency_ff5 = tan_ff5,
  ff5_factor_stats = factor_points$ff5,
  xlim = mv_limits$xlim,
  ylim = mv_limits$ylim
)

### 7.18K Plot 2: Six-Factor ###################################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot_mv_ff6(
  frontier_ff6 = frontier_ff6,
  tangency_ff6 = tan_ff6,
  ff6_factor_stats = factor_points$ff6,
  xlim = mv_limits$xlim,
  ylim = mv_limits$ylim
)

### 7.18L Plot 3: Combined #####################################################
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot_mv_both(
  frontier_ff5 = frontier_ff5,
  frontier_ff6 = frontier_ff6,
  tangency_ff5 = tan_ff5,
  tangency_ff6 = tan_ff6,
  ff5_factor_stats = factor_points$ff5,
  ff6_factor_stats = factor_points$ff6,
  xlim = mv_limits$xlim,
  ylim = mv_limits$ylim
)

### 7.18M Plot Factor Loadings #################################################
### Plot FF5 factor loadings across chi deciles ################################

library(data.table)

### 1. Load the saved FF5 loading table ########################################
ff5_loading_tables <- readRDS("table_ff5_alpha_loadings_chi_sic_baseline.rds")

tab <- ff5_loading_tables[["CHI_SIC_DEC10_BASELINE"]]

### 2. Keep only decile columns (not SMQ) ######################################
decile_cols <- paste0("P", 1:10)
factor_rows <- c("MKT", "SMB", "HML", "RMW", "CMA")

tab_dt <- as.data.table(tab, keep.rownames = "Factor")
tab_dt <- tab_dt[Factor %in% factor_rows, c("Factor", decile_cols), with = FALSE]

### 3. Extract numeric estimates from cells like "0.95***\n(4.21)" ############
extract_estimate <- function(x) {
  x <- sub("\n.*", "", x)
  x <- gsub("\\*", "", x)
  as.numeric(x)
}

plot_dt <- melt(
  tab_dt,
  id.vars = "Factor",
  variable.name = "Decile",
  value.name = "Cell"
)

plot_dt[, Loading := vapply(Cell, extract_estimate, numeric(1))]
plot_dt[, Decile_num := as.integer(gsub("P", "", Decile))]

setorder(plot_dt, Factor, Decile_num)


### 7.18N Plot FF5 factor loadings across chi deciles ##########################

par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

plot(
  x = 1:10,
  y = plot_dt[Factor == "MKT", Loading],
  type = "b",
  pch = 16,
  lwd = 2,
  ylim = range(plot_dt$Loading, na.rm = TRUE),
  xlab = "Chi decile portfolio",
  ylab = "Factor loading",
  main = "FF5 Factor Loadings Across Chi-Sorted Deciles",
  xaxt = "n"
)

axis(1, at = 1:10, labels = paste0("P", 1:10))
abline(h = 0, lty = 2)

lines(1:10, plot_dt[Factor == "SMB", Loading], type = "b", pch = 16, lwd = 2, col = 2)
lines(1:10, plot_dt[Factor == "HML", Loading], type = "b", pch = 16, lwd = 2, col = 3)
lines(1:10, plot_dt[Factor == "RMW", Loading], type = "b", pch = 16, lwd = 2, col = 4)
lines(1:10, plot_dt[Factor == "CMA", Loading], type = "b", pch = 16, lwd = 2, col = 5)

legend(
  "topright",
  legend = c("MKT", "SMB", "HML", "RMW", "CMA"),
  col = 1:5,
  lty = 1,
  pch = 16,
  lwd = 2,
  bty = "n"
)


### 7.18N Plot FF5 factor loadings and alpha as bar charts ####################
### 1. Include alpha as well
factor_rows <- c("FF5_alpha", "MKT", "SMB", "HML", "RMW", "CMA")

### 2. Rebuild plotting table with alpha included
tab_dt <- as.data.table(tab, keep.rownames = "Factor")
tab_dt <- tab_dt[Factor %in% factor_rows, c("Factor", decile_cols), with = FALSE]

plot_dt <- melt(
  tab_dt,
  id.vars = "Factor",
  variable.name = "Decile",
  value.name = "Cell"
)

plot_dt[, Loading := vapply(Cell, extract_estimate, numeric(1))]
plot_dt[, Decile_num := as.integer(gsub("P", "", Decile))]

setorder(plot_dt, Factor, Decile_num)

### 3. Plot
par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))

for (fac in factor_rows) {
  y <- plot_dt[Factor == fac, Loading]

  bp <- barplot(
    height = y,
    names.arg = paste0("P", 1:10),
    col = "#08306B",
    border = "black",
    main = fac,
    xlab = "Chi decile portfolio",
    ylab = ifelse(fac == "FF5_alpha", "Alpha (% per month)", "Loading"),
    las = 1
  )

  abline(h = 0, lty = 2)
}

### 7.18N Plot FF5 factor loadings and alpha as bar charts ####################

pdf(
  file   = file.path(FIGURES_DIR, "figure_ff5_alpha_and_loadings_barplots.pdf"),
  width  = 10,
  height = 8
)

### 1. Include alpha as well
factor_rows <- c("FF5_alpha", "MKT", "SMB", "HML", "RMW", "CMA")

### 2. Rebuild plotting table with alpha included
tab_dt <- as.data.table(tab, keep.rownames = "Factor")
tab_dt <- tab_dt[Factor %in% factor_rows, c("Factor", decile_cols), with = FALSE]

plot_dt <- melt(
  tab_dt,
  id.vars = "Factor",
  variable.name = "Decile",
  value.name = "Cell"
)

plot_dt[, Loading := vapply(Cell, extract_estimate, numeric(1))]
plot_dt[, Decile_num := as.integer(gsub("P", "", Decile))]

setorder(plot_dt, Factor, Decile_num)

### 3. Plot
par(mfrow = c(3, 2), mar = c(4, 4, 2.5, 1))

for (fac in factor_rows) {
  y <- plot_dt[Factor == fac, Loading]

  barplot(
    height = y,
    names.arg = paste0("P", 1:10),
    col = "#08306B",
    border = "black",
    main = fac,
    xlab = "Chi decile portfolio",
    ylab = ifelse(fac == "FF5_alpha", "Alpha (% per month)", "Loading"),
    las = 1
  )

  abline(h = 0, lty = 2)
}

dev.off()

### 7.18N Plot FF5 factor loadings as grouped bar chart ########################

loading_mat <- dcast(
  plot_dt,
  Factor ~ Decile,
  value.var = "Loading"
)

loading_mat_num <- as.matrix(loading_mat[, -1])
rownames(loading_mat_num) <- loading_mat$Factor

factor_cols <- c(
  "gray40", # MKT
  "steelblue3", # SMB
  "purple3", # HML
  "mediumpurple4", # RMW
  "orchid4" # CMA
)

par(mfrow = c(1, 1), mar = c(5, 4, 3, 1))

barplot(
  loading_mat_num,
  beside = TRUE,
  col = factor_cols,
  names.arg = paste0("P", 1:10),
  xlab = "Chi decile portfolio",
  ylab = "Factor loading",
  main = "FF5 Factor Loadings Across Chi-Sorted Deciles",
  legend.text = rownames(loading_mat_num),
  args.legend = list(x = "topright", bty = "n", fill = factor_cols)
)

abline(h = 0, lty = 2)
