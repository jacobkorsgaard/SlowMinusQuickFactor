## 10. Robustness Tests #######################################################
cat("\014")
rm(list = ls())
graphics.off()

### 10.1 Setup ################################################################
library(data.table)
library(lubridate)
library(lmtest)
library(sandwich)
library(zoo)

source(file.path("R", "00_config.R"))

START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")
NW_LAGS <- 12
ROLL_WIN <- 60

### 10.2 Load Data ############################################################
panel <- readRDS("master_factors_testassets.rds")
setDT(panel)
panel[, mdate := as.Date(mdate)]

factor_test_panel <- copy(panel)

panel <- panel[
  mdate >= START_DATE &
    mdate <= END_DATE
]

factor_test_panel <- factor_test_panel[
  mdate >= START_DATE &
    mdate <= END_DATE
]

raw_data <- readRDS("crsp_compustat_with_adj_cost.rds")
setDT(raw_data)
raw_data[, mdate := as.Date(mdate)]

raw_data <- raw_data[
  mdate >= START_DATE &
    mdate <= END_DATE
]

raw_data[, year := year(mdate)]
raw_data[, month := month(mdate)]
raw_data[, ffyear := fifelse(month >= 7, year, year - 1)]

if (!"sic2" %in% names(raw_data)) {
  raw_data[, sic2 := fifelse(!is.na(sic), as.integer(sic) %/% 100L, NA_integer_)]
}

ff5 <- readRDS("ff5_factors_internal.rds")
setDT(ff5)
ff5[, mdate := as.Date(mdate)]
ff5 <- ff5[
  mdate >= START_DATE &
    mdate <= END_DATE
]

SMQ <- readRDS("SMQ_factors.rds")
setDT(SMQ)
SMQ[, mdate := as.Date(mdate)]
SMQ <- SMQ[
  mdate >= START_DATE &
    mdate <= END_DATE
]

labelled_data <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(labelled_data)
labelled_data[, mdate := as.Date(mdate)]
labelled_data <- labelled_data[
  mdate >= START_DATE &
    mdate <= END_DATE
]

rf_monthly <- unique(factor_test_panel[, .(mdate, RF)])

### 10.3 Helpers ##############################################################
star_p <- function(p) {
  fifelse(
    is.na(p), "",
    fifelse(
      p < 0.01, "***",
      fifelse(
        p < 0.05, "**",
        fifelse(p < 0.10, "*", "")
      )
    )
  )
}

star_t <- function(t) {
  if (is.na(t)) "" else if (abs(t) >= 2.58) "***" else if (abs(t) >= 1.96) "**" else if (abs(t) >= 1.65) "*" else ""
}

round_dt <- function(dt, digits = 4) {
  out <- copy(dt)
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  out
}

cum_safe <- function(x) {
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- cumprod(1 + x[ok]) - 1
  out
}

ff_30_70_breaks <- function(x) {
  qs <- quantile(x, probs = c(0.3, 0.7), na.rm = TRUE, type = 1)
  c(low = as.numeric(qs[1]), high = as.numeric(qs[2]))
}

nw_tstat_mean <- function(x, lags = NW_LAGS) {
  x <- x[!is.na(x)]
  if (length(x) < 24) {
    return(NA_real_)
  }
  fit <- lm(x ~ 1)
  nw <- NeweyWest(fit, lag = min(lags, floor(length(x) / 4)), prewhite = FALSE, adjust = TRUE)
  ct <- coeftest(fit, vcov. = nw)
  as.numeric(ct[1, "t value"])
}

mean_nw <- function(y, lags = NW_LAGS) {
  dt <- data.table(y = y)
  dt <- dt[complete.cases(dt)]

  if (nrow(dt) < 24) {
    return(data.table(
      estimate = NA_real_,
      tstat    = NA_real_,
      pval     = NA_real_,
      stars    = ""
    ))
  }

  fit <- lm(y ~ 1, data = dt)
  nw <- NeweyWest(
    fit,
    lag = min(lags, floor(nrow(dt) / 4)),
    prewhite = FALSE,
    adjust = TRUE
  )
  ct <- coeftest(fit, vcov. = nw)

  data.table(
    estimate = as.numeric(ct[1, 1]),
    tstat    = as.numeric(ct[1, 3]),
    pval     = as.numeric(ct[1, 4]),
    stars    = star_p(as.numeric(ct[1, 4]))
  )
}

nw_mean_test <- function(x, lags = NW_LAGS) {
  tmp <- data.table(x = x)
  tmp <- tmp[complete.cases(tmp)]

  if (nrow(tmp) < 24) {
    return(data.table(
      estimate = NA_real_,
      tstat    = NA_real_,
      pval     = NA_real_,
      stars    = ""
    ))
  }

  fit <- lm(x ~ 1, data = tmp)
  nw <- NeweyWest(
    fit,
    lag = min(lags, floor(nrow(tmp) / 4)),
    prewhite = FALSE,
    adjust = TRUE
  )
  ct <- coeftest(fit, vcov. = nw)

  data.table(
    estimate = as.numeric(ct[1, 1]),
    tstat    = as.numeric(ct[1, 3]),
    pval     = as.numeric(ct[1, 4]),
    stars    = star_p(as.numeric(ct[1, 4]))
  )
}

run_ts_reg_nw <- function(y, X, lags = NW_LAGS) {
  dt <- data.table(y = y)
  dt <- cbind(dt, as.data.table(X))
  dt <- dt[complete.cases(dt)]

  if (nrow(dt) < 24) {
    return(NULL)
  }

  fml <- as.formula(paste("y ~", paste(names(X), collapse = " + ")))
  fit <- lm(fml, data = dt)
  nw <- NeweyWest(
    fit,
    lag = min(lags, floor(nrow(dt) / 4)),
    prewhite = FALSE,
    adjust = TRUE
  )
  ct <- coeftest(fit, vcov. = nw)
  fit_sum <- summary(fit)

  out <- data.table(
    term     = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    tstat    = as.numeric(ct[, 3]),
    pval     = as.numeric(ct[, 4])
  )
  out[, stars := star_p(pval)]
  out[, adj_r2 := fit_sum$adj.r.squared]
  out[]
}

run_ts_reg_nw_formula <- function(formula, dt, lags = NW_LAGS) {
  tmp <- copy(dt)
  vars_needed <- all.vars(formula)
  tmp <- tmp[, ..vars_needed]
  tmp <- tmp[complete.cases(tmp)]

  if (nrow(tmp) < 24) {
    return(NULL)
  }

  fit <- lm(formula, data = tmp)
  nw <- NeweyWest(
    fit,
    lag = min(lags, floor(nrow(tmp) / 4)),
    prewhite = FALSE,
    adjust = TRUE
  )
  ct <- coeftest(fit, vcov. = nw)
  fit_sum <- summary(fit)

  out <- data.table(
    term     = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    tstat    = as.numeric(ct[, 3]),
    pval     = as.numeric(ct[, 4]),
    stars    = star_p(as.numeric(ct[, 4]))
  )
  out[, adj_r2 := fit_sum$adj.r.squared]
  out
}

format_cell <- function(est, tstat, stars, digits_est = 2, digits_t = 2, scale = 100) {
  est_s <- scale * est
  sprintf(
    paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
    est_s, stars, tstat
  )
}

format_cell_with_p <- function(est, tstat, pval, stars, digits_est = 2, digits_t = 2, digits_p = 3, scale = 100) {
  est_s <- scale * est
  sprintf(
    paste0("%.", digits_est, "f%s\n(t=%.", digits_t, "f, p=%.", digits_p, "f)"),
    est_s, stars, tstat, pval
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

roll_mean_safe <- function(x, k = 60) {
  zoo::rollapply(
    x,
    width = k,
    FUN = mean,
    fill = NA,
    align = "right",
    na.rm = TRUE
  )
}

roll_sharpe_safe <- function(x, k = 60) {
  zoo::rollapply(
    x,
    width = k,
    FUN = function(z) {
      z <- z[!is.na(z)]
      if (length(z) < 24) {
        return(NA_real_)
      }
      s <- sd(z)
      if (is.na(s) || s == 0) {
        return(NA_real_)
      }
      mean(z) / s * sqrt(12)
    },
    fill = NA,
    align = "right"
  )
}

compute_hhi <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  shares <- as.numeric(table(x)) / length(x)
  sum(shares^2)
}

compute_weighted_hhi <- function(industry, weight) {
  dt <- data.table(industry = industry, weight = weight)
  dt <- dt[
    !is.na(industry) &
      !is.na(weight) &
      is.finite(weight) &
      weight > 0
  ]

  if (nrow(dt) == 0) {
    return(NA_real_)
  }

  ind_w <- dt[, .(industry_weight = sum(weight, na.rm = TRUE)), by = industry]
  total_w <- ind_w[, sum(industry_weight, na.rm = TRUE)]

  if (is.na(total_w) || total_w <= 0) {
    return(NA_real_)
  }

  ind_w[, share := industry_weight / total_w]
  sum(ind_w$share^2)
}

effective_n <- function(hhi) {
  ifelse(is.na(hhi) | hhi <= 0, NA_real_, 1 / hhi)
}

make_leg_hhi_table <- function(dt, grp_var, factor_name, industry_var = "sic2") {
  tmp <- dt[
    !is.na(get(grp_var)) &
      !is.na(get(industry_var)),
    .(
      N_firms = uniqueN(permno),
      N_industries = uniqueN(get(industry_var)),
      HHI = compute_hhi(get(industry_var))
    ),
    by = .(ffyear, Leg = get(grp_var))
  ]

  tmp[, Effective_N := effective_n(HHI)]
  tmp[, Factor := factor_name]
  setcolorder(tmp, c("Factor", "ffyear", "Leg", "N_firms", "N_industries", "HHI", "Effective_N"))
  tmp[]
}

make_leg_weighted_hhi_table <- function(dt, grp_var, factor_name, industry_var = "sic2", weight_var = "me_june") {
  tmp <- dt[
    !is.na(get(grp_var)) &
      !is.na(get(industry_var)) &
      !is.na(get(weight_var)) &
      is.finite(get(weight_var)) &
      get(weight_var) > 0,
    .(
      N_firms = uniqueN(permno),
      N_industries = uniqueN(get(industry_var)),
      Total_weight = sum(get(weight_var), na.rm = TRUE),
      HHI_weighted = compute_weighted_hhi(get(industry_var), get(weight_var))
    ),
    by = .(ffyear, Leg = get(grp_var))
  ]

  tmp[, Effective_N_weighted := effective_n(HHI_weighted)]
  tmp[, Factor := factor_name]
  setcolorder(
    tmp,
    c("Factor", "ffyear", "Leg", "N_firms", "N_industries", "Total_weight", "HHI_weighted", "Effective_N_weighted")
  )
  tmp[]
}

make_leg_hhi_summary <- function(hhi_dt, hhi_col = "HHI", effn_col = "Effective_N") {
  hhi_dt[
    ,
    .(
      Avg_HHI          = mean(get(hhi_col), na.rm = TRUE),
      Median_HHI       = median(get(hhi_col), na.rm = TRUE),
      Min_HHI          = min(get(hhi_col), na.rm = TRUE),
      Max_HHI          = max(get(hhi_col), na.rm = TRUE),
      Avg_Effective_N  = mean(get(effn_col), na.rm = TRUE),
      Avg_N_firms      = mean(N_firms, na.rm = TRUE),
      Avg_N_industries = mean(N_industries, na.rm = TRUE)
    ),
    by = .(Factor, Leg)
  ][order(Factor, Leg)]
}

make_summary_table <- function(dt, factor_vec) {
  means <- sapply(factor_vec, function(f) 100 * mean(dt[[f]], na.rm = TRUE))
  sds <- sapply(factor_vec, function(f) 100 * sd(dt[[f]], na.rm = TRUE))
  sharpe <- sapply(factor_vec, function(f) {
    m <- mean(dt[[f]], na.rm = TRUE)
    s <- sd(dt[[f]], na.rm = TRUE)
    if (is.na(s) || s == 0) {
      return(NA_real_)
    }
    (m / s) * sqrt(12)
  })
  tstats <- sapply(factor_vec, function(f) nw_tstat_mean(dt[[f]], NW_LAGS))
  stars <- sapply(tstats, star_t)

  tab <- data.table(
    Statistic = c("Mean", "Std. Dev.", "Sharpe Ratio", "t-Statistic")
  )

  for (f in factor_vec) {
    tab[, (f) := c(
      sprintf("%.2f%s", means[f], stars[f]),
      sprintf("%.2f", sds[f]),
      sprintf("%.2f", sharpe[f]),
      sprintf("%.2f", tstats[f])
    )]
  }

  tab
}

make_factor_table <- function(dt, fac_names, lags = NW_LAGS, digits_est = 2, digits_t = 2, scale_est = 100) {
  res <- list()

  for (y in fac_names) {
    x <- setdiff(fac_names, y)
    reg <- run_ts_reg_nw_formula(as.formula(paste(y, "~", paste(x, collapse = " + "))), dt, lags)
    if (is.null(reg)) next

    row <- data.table(dep = y)
    a <- reg[term == "(Intercept)"]

    row[, Int := sprintf(
      paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
      scale_est * a$estimate, a$stars, a$tstat
    )]

    for (xx in x) {
      b <- reg[term == xx]
      row[, (xx) := sprintf(
        paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
        scale_est * b$estimate, b$stars, b$tstat
      )]
    }

    row[, Adj_R2 := sprintf("%.4f", unique(reg$adj_r2))]
    res[[y]] <- row
  }

  tab <- rbindlist(res, fill = TRUE)
  col_order <- c("dep", "Int", fac_names, "Adj_R2")
  col_order <- col_order[col_order %in% names(tab)]
  tab[, ..col_order]
}

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

  list(
    stats = out[, .(Factor, Leg, Mean_excess, Tstat, Vol, Sharpe)],
    spread_series = spread[, .(mdate, exret)]
  )
}

display_all_5x5 <- function(mat_list, title, digits = 3) {
  cat("\n=====================================\n")
  cat(title, "\n")
  cat("=====================================\n\n")
  for (name in names(mat_list)) {
    cat("-----", name, "-----\n")
    print(round(mat_list[[name]], digits))
    cat("\n")
  }
}

display_2x4x4_tables <- function(tab_list, title, digits = 3) {
  cat("\n============================================================\n")
  cat(title, "\n")
  cat("============================================================\n\n")

  for (nm in names(tab_list)) {
    cat("------------------------------------------------------------\n")
    cat(nm, "\n")
    cat("------------------------------------------------------------\n")
    cat("\nSmall\n")
    print(round(tab_list[[nm]]$S, digits))
    cat("\nBig\n")
    print(round(tab_list[[nm]]$B, digits))
    cat("\n")
  }
}

make_subperiod_summary_table <- function(dt, factor_vec, lags = NW_LAGS) {
  out <- data.table(
    Statistic = c("Mean", "t-Statistic", "Std. Dev.", "Sharpe Ratio")
  )

  for (f in factor_vec) {
    x <- dt[[f]]
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    t <- nw_tstat_mean(x, lags)
    sr <- ifelse(is.na(s) || s == 0, NA_real_, m / s * sqrt(12))
    stars <- star_t(t)

    out[, (f) := c(
      sprintf("%.2f%s", 100 * m, stars),
      sprintf("%.2f", t),
      sprintf("%.2f", 100 * s),
      sprintf("%.2f", sr)
    )]
  }

  out
}

make_subperiod_factor_regression_table <- function(dt, smq_col, lags = NW_LAGS, digits_est = 2, digits_t = 2, scale_est = 100) {
  stopifnot(all(c("HML", "CMA", smq_col) %in% names(dt)))

  facs <- c("HML", "CMA", smq_col)

  reg_hml <- run_ts_reg_nw_formula(
    as.formula(paste("HML ~", paste(setdiff(facs, "HML"), collapse = " + "))),
    dt,
    lags
  )

  reg_cma <- run_ts_reg_nw_formula(
    as.formula(paste("CMA ~", paste(setdiff(facs, "CMA"), collapse = " + "))),
    dt,
    lags
  )

  reg_smq <- run_ts_reg_nw_formula(
    as.formula(paste(smq_col, "~", paste(setdiff(facs, smq_col), collapse = " + "))),
    dt,
    lags
  )

  build_panel <- function(reg, dep_name, regressor_order) {
    row_int <- reg[term == "(Intercept)"]

    out <- data.table(
      term = c("Int.", regressor_order, "Adj. R2")
    )

    out[, value := ""]

    out[term == "Int.", value := sprintf(
      paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
      scale_est * row_int$estimate,
      row_int$stars,
      row_int$tstat
    )]

    for (rr in regressor_order) {
      row_rr <- reg[term == rr]
      if (nrow(row_rr) > 0) {
        out[term == rr, value := sprintf(
          paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
          scale_est * row_rr$estimate,
          row_rr$stars,
          row_rr$tstat
        )]
      }
    }

    out[term == "Adj. R2", value := sprintf("%.4f", unique(reg$adj_r2))]
    out[, dep := dep_name]
    out
  }

  panel_hml <- build_panel(reg_hml, "HML", c("CMA", smq_col))
  panel_cma <- build_panel(reg_cma, "CMA", c("HML", smq_col))
  panel_smq <- build_panel(reg_smq, "SMQ", c("HML", "CMA"))

  list(
    HML = panel_hml,
    CMA = panel_cma,
    SMQ = panel_smq
  )
}

run_alpha_model_on_series <- function(dt, dep_col, rhs_cols, lags = NW_LAGS) {
  stopifnot(dep_col %in% names(dt))
  stopifnot(all(rhs_cols %in% names(dt)))

  reg <- run_ts_reg_nw(
    y = dt[[dep_col]],
    X = dt[, ..rhs_cols],
    lags = lags
  )

  if (is.null(reg)) {
    return(data.table(
      model = paste(rhs_cols, collapse = "+"),
      alpha = NA_real_,
      alpha_t = NA_real_,
      alpha_p = NA_real_,
      alpha_stars = "",
      adj_r2 = NA_real_
    ))
  }

  a <- reg[term == "(Intercept)"][1]

  data.table(
    model = paste(rhs_cols, collapse = "+"),
    alpha = a$estimate,
    alpha_t = a$tstat,
    alpha_p = a$pval,
    alpha_stars = a$stars,
    adj_r2 = a$adj_r2
  )
}

make_smq_comparison_table <- function(dt, smq_cols, lags = NW_LAGS) {
  out <- list()

  for (smq_col in smq_cols) {
    stopifnot(smq_col %in% names(dt))

    mean_test <- nw_mean_test(dt[[smq_col]], lags)
    vol_monthly <- sd(dt[[smq_col]], na.rm = TRUE)
    sharpe <- ifelse(
      is.na(sd(dt[[smq_col]], na.rm = TRUE)) || sd(dt[[smq_col]], na.rm = TRUE) == 0,
      NA_real_,
      mean(dt[[smq_col]], na.rm = TRUE) / sd(dt[[smq_col]], na.rm = TRUE) * sqrt(12)
    )

    capm <- run_alpha_model_on_series(dt, smq_col, c("MKT"), lags)
    ff3m <- run_alpha_model_on_series(dt, smq_col, c("MKT", "SMB", "HML"), lags)
    ff5m <- run_alpha_model_on_series(dt, smq_col, c("MKT", "SMB", "HML", "RMW", "CMA"), lags)

    row <- data.table(
      factor = smq_col,
      mean_excess = mean_test$estimate,
      mean_tstat = mean_test$tstat,
      mean_pval = mean_test$pval,
      mean_stars = mean_test$stars,
      vol_monthly = vol_monthly,
      sharpe = sharpe,
      CAPM_alpha = capm$alpha,
      CAPM_tstat = capm$alpha_t,
      CAPM_pval = capm$alpha_p,
      CAPM_stars = capm$alpha_stars,
      CAPM_adj_r2 = capm$adj_r2,
      FF3_alpha = ff3m$alpha,
      FF3_tstat = ff3m$alpha_t,
      FF3_pval = ff3m$alpha_p,
      FF3_stars = ff3m$alpha_stars,
      FF3_adj_r2 = ff3m$adj_r2,
      FF5_alpha = ff5m$alpha,
      FF5_tstat = ff5m$alpha_t,
      FF5_pval = ff5m$alpha_p,
      FF5_stars = ff5m$alpha_stars,
      FF5_adj_r2 = ff5m$adj_r2
    )

    out[[smq_col]] <- row
  }

  rbindlist(out, fill = TRUE)
}

make_smq_display_table <- function(dt) {
  out <- copy(dt)

  out[, Mean := sprintf(
    "%.2f%s\n(%.2f, p=%.3f)",
    100 * mean_excess, mean_stars, mean_tstat, mean_pval
  )]
  out[, Volatility := sprintf("%.2f", 100 * vol_monthly)]
  out[, Sharpe := sprintf("%.2f", sharpe)]

  out[, CAPM := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * CAPM_alpha, CAPM_stars, CAPM_tstat, CAPM_pval
  )]
  out[, FF3 := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * FF3_alpha, FF3_stars, FF3_tstat, FF3_pval
  )]
  out[, FF5 := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * FF5_alpha, FF5_stars, FF5_tstat, FF5_pval
  )]

  out[, .(factor, Mean, Volatility, Sharpe, CAPM, FF3, FF5)]
}

extract_alpha_table <- function(reg_dt, portfolio_cols = c("portfolio", "size", "char", "i", "j")) {
  stopifnot("term" %in% names(reg_dt))
  out <- copy(reg_dt[term == "(Intercept)"])
  keep <- unique(c(portfolio_cols, "estimate", "tstat", "pval", "stars", "adj_r2"))
  keep <- keep[keep %in% names(out)]
  out <- out[, ..keep]
  setnames(out,
    old = c("estimate", "tstat", "pval", "stars"),
    new = c("alpha", "alpha_tstat", "alpha_pval", "alpha_stars")
  )
  out[]
}

summarize_significant_alphas <- function(alpha_dt, family_name, chi_pos = NULL, size_pos = NULL, levels = c(0.10, 0.05, 0.01)) {
  out <- list()

  base_row <- data.table(
    family = family_name,
    scope = "overall",
    bucket = "all",
    N = nrow(alpha_dt)
  )

  for (lvl in levels) {
    sig <- alpha_dt$alpha_pval < lvl
    base_row[, paste0("n_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig, na.rm = TRUE)]
    base_row[, paste0("prop_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := mean(sig, na.rm = TRUE)]
    base_row[, paste0("n_pos_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & alpha_dt$alpha > 0, na.rm = TRUE)]
    base_row[, paste0("n_neg_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & alpha_dt$alpha < 0, na.rm = TRUE)]
  }
  out[["overall"]] <- base_row

  if (!is.null(size_pos) && size_pos %in% names(alpha_dt)) {
    out_size <- rbindlist(lapply(sort(unique(alpha_dt[[size_pos]])), function(ss) {
      tmp <- alpha_dt[get(size_pos) == ss]
      row <- data.table(
        family = family_name,
        scope = "size",
        bucket = as.character(ss),
        N = nrow(tmp)
      )
      for (lvl in levels) {
        sig <- tmp$alpha_pval < lvl
        row[, paste0("n_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig, na.rm = TRUE)]
        row[, paste0("prop_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := mean(sig, na.rm = TRUE)]
        row[, paste0("n_pos_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & tmp$alpha > 0, na.rm = TRUE)]
        row[, paste0("n_neg_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & tmp$alpha < 0, na.rm = TRUE)]
      }
      row
    }), fill = TRUE)
    out[["size"]] <- out_size
  }

  if (!is.null(chi_pos) && chi_pos %in% names(alpha_dt)) {
    out_chi <- rbindlist(lapply(sort(unique(alpha_dt[[chi_pos]])), function(cc) {
      tmp <- alpha_dt[get(chi_pos) == cc]
      row <- data.table(
        family = family_name,
        scope = "chi_bucket",
        bucket = as.character(cc),
        N = nrow(tmp)
      )
      for (lvl in levels) {
        sig <- tmp$alpha_pval < lvl
        row[, paste0("n_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig, na.rm = TRUE)]
        row[, paste0("prop_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := mean(sig, na.rm = TRUE)]
        row[, paste0("n_pos_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & tmp$alpha > 0, na.rm = TRUE)]
        row[, paste0("n_neg_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & tmp$alpha < 0, na.rm = TRUE)]
      }
      row
    }), fill = TRUE)
    out[["chi"]] <- out_chi
  }

  if (!is.null(size_pos) && !is.null(chi_pos) && size_pos %in% names(alpha_dt) && chi_pos %in% names(alpha_dt)) {
    out_size_chi <- alpha_dt[
      ,
      {
        row <- data.table(N = .N)
        for (lvl in levels) {
          sig <- alpha_pval < lvl
          row[, paste0("n_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig, na.rm = TRUE)]
          row[, paste0("prop_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := mean(sig, na.rm = TRUE)]
          row[, paste0("n_pos_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & alpha > 0, na.rm = TRUE)]
          row[, paste0("n_neg_sig_", sub("\\.", "", sprintf("%.2f", lvl))) := sum(sig & alpha < 0, na.rm = TRUE)]
        }
        row
      },
      by = c(size_pos, chi_pos)
    ]
    out_size_chi[, family := family_name]
    out_size_chi[, scope := "size_x_chi"]
    setnames(out_size_chi, c(size_pos, chi_pos), c("size_bucket", "chi_bucket"))
    out[["size_x_chi"]] <- out_size_chi
  }

  rbindlist(out, fill = TRUE)
}

make_alpha_count_comparison <- function(alpha_list_named, family_name, size_col = NULL, chi_col = NULL) {
  out <- list()

  for (model_name in names(alpha_list_named)) {
    alpha_dt <- copy(alpha_list_named[[model_name]])
    alpha_dt[, model := model_name]

    tmp <- summarize_significant_alphas(
      alpha_dt = alpha_dt,
      family_name = family_name,
      chi_pos = chi_col,
      size_pos = size_col
    )
    tmp[, model := model_name]
    out[[model_name]] <- tmp
  }

  rbindlist(out, fill = TRUE)
}

wide_5x5_to_long_regex <- function(panel, regex_pattern, strip_prefix) {
  cols <- grep(regex_pattern, names(panel), value = TRUE)
  stopifnot(length(cols) > 0)

  long <- melt(
    panel[, c("mdate", cols), with = FALSE],
    id.vars = "mdate",
    variable.name = "portfolio",
    value.name = "vwret"
  )
  long <- long[!is.na(vwret)]

  parse_5x5 <- function(col, strip_prefix) {
    x <- sub(paste0("^", strip_prefix), "", col)
    m <- regexec("^(S[1-5])([A-Z])(\\d)$", x)
    g <- regmatches(x, m)[[1]]
    if (length(g) != 4) {
      return(NULL)
    }
    data.table(size5 = g[2], char5 = paste0(g[3], g[4]))
  }

  parsed <- rbindlist(lapply(long$portfolio, parse_5x5, strip_prefix = strip_prefix), fill = TRUE)
  cbind(long, parsed)
}

make_5x5_mean_matrix <- function(port_dt, size_var = "size5", char_var, ret_var = "vwret") {
  means <- port_dt[, .(mean_ret = mean(get(ret_var), na.rm = TRUE)), by = c(size_var, char_var)]
  means[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  means[, (char_var) := factor(get(char_var), levels = sort(unique(get(char_var))))]

  mat <- dcast(means, as.formula(paste(size_var, "~", char_var)), value.var = "mean_ret")
  mat_out <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat_out) <- mat[[size_var]]
  100 * mat_out
}

make_5x5_sd_matrix <- function(port_dt, size_var = "size5", char_var, ret_var = "vwret_excess") {
  sds <- port_dt[, .(sd_ret = sd(get(ret_var), na.rm = TRUE)), by = c(size_var, char_var)]
  sds[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  sds[, (char_var) := factor(get(char_var), levels = sort(unique(get(char_var))))]

  mat <- dcast(sds, as.formula(paste(size_var, "~", char_var)), value.var = "sd_ret")
  mat_out <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat_out) <- mat[[size_var]]
  100 * mat_out
}

make_5x5_sr_matrix <- function(port_dt, size_var = "size5", char_var, ret_var = "vwret_excess") {
  stats <- port_dt[, .(
    mean_ret = mean(get(ret_var), na.rm = TRUE),
    sd_ret   = sd(get(ret_var), na.rm = TRUE)
  ), by = c(size_var, char_var)]

  stats[, sr := mean_ret / sd_ret]
  stats[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  stats[, (char_var) := factor(get(char_var), levels = sort(unique(get(char_var))))]

  mat <- dcast(stats, as.formula(paste(size_var, "~", char_var)), value.var = "sr")
  mat_out <- as.matrix(mat[, -1, with = FALSE])
  rownames(mat_out) <- mat[[size_var]]
  mat_out
}

make_excess <- function(port_dt, rf_dt) {
  tmp <- merge(port_dt, rf_dt, by = "mdate", all.x = TRUE)
  tmp[, vwret_excess := vwret - RF]
  tmp
}

filter_window <- function(dt, start_date = START_DATE, end_date = END_DATE) {
  dt[mdate >= start_date & mdate <= end_date]
}

build_excess_stats <- function(port_dt, char_var, rf_dt) {
  ex <- make_excess(filter_window(port_dt), rf_dt)
  list(
    mean = make_5x5_mean_matrix(ex, "size5", char_var, "vwret_excess"),
    sd   = make_5x5_sd_matrix(ex, "size5", char_var, "vwret_excess"),
    sr   = make_5x5_sr_matrix(ex, "size5", char_var, "vwret_excess")
  )
}

parse_5x5_from_col <- function(col, strip_prefix) {
  x <- sub(paste0("^", strip_prefix), "", col)
  m <- regexec("^(S[1-5])([A-Z])(\\d)$", x)
  g <- regmatches(x, m)[[1]]
  if (length(g) != 4) {
    return(NULL)
  }
  data.table(size = g[2], char = paste0(g[3], g[4]))
}

term_to_5x5_matrix <- function(reg_dt, term_name, digits_est = 2, digits_t = 2, scale_est = 100) {
  tmp <- reg_dt[term == term_name, .(size, char, estimate, tstat, stars)]
  if (nrow(tmp) == 0) {
    return(matrix("", 5, 5))
  }

  tmp[, est_s := scale_est * estimate]
  tmp[, cell := sprintf(
    paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
    est_s, stars, tstat
  )]

  tmp[, size := factor(size, levels = paste0("S", 1:5))]
  tmp[, char := factor(char, levels = sort(unique(char)))]

  mat <- dcast(tmp, size ~ char, value.var = "cell", fill = "")
  out <- as.matrix(mat[, -1, with = FALSE])
  rownames(out) <- mat$size
  out
}

run_asset_block_regs_regex <- function(panel, regex_pattern, strip_prefix, rhs, lags = NW_LAGS) {
  asset_cols <- grep(regex_pattern, names(panel), value = TRUE)
  stopifnot(length(asset_cols) > 0)

  X <- panel[, ..rhs]
  res <- list()

  for (col in asset_cols) {
    lab <- parse_5x5_from_col(col, strip_prefix)
    if (is.null(lab)) next

    y <- panel[[col]] - panel[["RF"]]
    est <- run_ts_reg_nw(y, X, lags)
    if (is.null(est)) next

    est[, portfolio := col]
    est[, size := lab$size]
    est[, char := lab$char]
    res[[col]] <- est
  }

  rbindlist(res, fill = TRUE)
}

wide_2x4x4_to_long_regex <- function(panel, regex_pattern, strip_prefix, expected_n = 32) {
  cols <- grep(regex_pattern, names(panel), value = TRUE)

  if (length(cols) == 0) {
    return(NULL)
  }

  if (length(cols) != expected_n) {
    warning(
      sprintf(
        "Pattern %s matched %d columns, not %d",
        regex_pattern, length(cols), expected_n
      )
    )
    return(NULL)
  }

  long <- melt(
    data = panel[, c("mdate", cols), with = FALSE],
    id.vars = "mdate",
    variable.name = "portfolio",
    value.name = "vwret"
  )

  long <- long[!is.na(vwret)]

  parse_2x4x4 <- function(col, strip_prefix) {
    x <- sub(paste0("^", strip_prefix), "", col)
    m <- regexec("^([SB])(\\d)(\\d)$", x)
    g <- regmatches(x, m)[[1]]

    if (length(g) != 4) {
      return(data.table(size = NA_character_, i = NA_character_, j = NA_character_))
    }

    data.table(
      size = g[2],
      i = g[3],
      j = g[4]
    )
  }

  parsed <- rbindlist(
    lapply(long$portfolio, parse_2x4x4, strip_prefix = strip_prefix),
    fill = TRUE
  )

  long[, c("size", "i", "j") := parsed[, .(size, i, j)]]

  setDT(long)
  long[]
}

make_mean_table_2x4x4 <- function(long_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  means <- x[, .(mean_ret = mean(vwret, na.rm = TRUE)), by = .(size, i, j)]

  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  means <- merge(grid, means, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(means[size == sz], i ~ j, value.var = "mean_ret", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_mean_excess_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  means <- x[, .(mean_excess = mean(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  means <- merge(grid, means, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(means[size == sz], i ~ j, value.var = "mean_excess", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sd_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  sds <- x[, .(sd_ret = sd(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  sds <- merge(grid, sds, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(sds[size == sz], i ~ j, value.var = "sd_ret", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sr_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  stats <- x[, .(
    mean_ret = mean(ret_excess, na.rm = TRUE),
    sd_ret   = sd(ret_excess, na.rm = TRUE)
  ), by = .(size, i, j)]

  stats[, sr := mean_ret / sd_ret]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  stats <- merge(grid, stats, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(stats[size == sz], i ~ j, value.var = "sr", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- mat_out
  }
  out
}

make_mean_excess_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  means <- x[, .(mean_excess = mean(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  means <- merge(grid, means, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(means[size == sz], i ~ j, value.var = "mean_excess", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sd_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  sds <- x[, .(sd_ret = sd(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  sds <- merge(grid, sds, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(sds[size == sz], i ~ j, value.var = "sd_ret", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sr_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- as.data.table(copy(long_dt))
  stopifnot(all(c("mdate", "size", "i", "j", "vwret") %in% names(x)))

  x <- x[get("mdate") >= start_date & get("mdate") <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  stats <- x[, .(
    mean_ret = mean(ret_excess, na.rm = TRUE),
    sd_ret   = sd(ret_excess, na.rm = TRUE)
  ), by = .(size, i, j)]

  stats[, sr := mean_ret / sd_ret]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  stats <- merge(grid, stats, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(stats[size == sz], i ~ j, value.var = "sr", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- mat_out
  }
  out
}

make_mean_excess_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- long_dt[mdate >= start_date & mdate <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  means <- x[, .(mean_excess = mean(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  means <- merge(grid, means, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(means[size == sz], i ~ j, value.var = "mean_excess", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sd_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- long_dt[mdate >= start_date & mdate <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  sds <- x[, .(sd_ret = sd(ret_excess, na.rm = TRUE)), by = .(size, i, j)]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  sds <- merge(grid, sds, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(sds[size == sz], i ~ j, value.var = "sd_ret", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- 100 * mat_out
  }
  out
}

make_sr_table_2x4x4 <- function(long_dt, rf_dt, start_date = START_DATE, end_date = END_DATE) {
  x <- long_dt[mdate >= start_date & mdate <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := vwret - RF]

  stats <- x[, .(
    mean_ret = mean(ret_excess, na.rm = TRUE),
    sd_ret   = sd(ret_excess, na.rm = TRUE)
  ), by = .(size, i, j)]

  stats[, sr := mean_ret / sd_ret]
  grid <- CJ(size = c("S", "B"), i = as.character(1:4), j = as.character(1:4), unique = TRUE)
  stats <- merge(grid, stats, by = c("size", "i", "j"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(stats[size == sz], i ~ j, value.var = "sr", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    colnames(mat_out) <- c("Q1", "Q2", "Q3", "Q4")
    out[[sz]] <- mat_out
  }
  out
}

parse_2x4x4_from_col <- function(col, strip_prefix) {
  x <- sub(paste0("^", strip_prefix), "", col)
  m <- regexec("^([SB])(\\d)(\\d)$", x)
  g <- regmatches(x, m)[[1]]
  if (length(g) != 4) {
    return(NULL)
  }
  data.table(size = g[2], i = g[3], j = g[4])
}

term_to_4x4_matrix <- function(reg_dt, term_name, size_flag, digits_est = 2, digits_t = 2, scale_est = 100) {
  tmp <- reg_dt[
    term == term_name & size == size_flag,
    .(i, j, estimate, tstat, stars)
  ]
  if (nrow(tmp) == 0) {
    return(matrix("", 4, 4))
  }

  tmp[, est_s := scale_est * estimate]
  tmp[, cell := sprintf(
    paste0("%.", digits_est, "f%s\n(%.", digits_t, "f)"),
    est_s, stars, tstat
  )]

  tmp[, i := factor(i, levels = c("1", "2", "3", "4"))]
  tmp[, j := factor(j, levels = c("1", "2", "3", "4"))]

  mat <- dcast(tmp, i ~ j, value.var = "cell", fill = "")
  out <- as.matrix(mat[, -1, with = FALSE])
  rownames(out) <- paste0("Q", mat$i)
  out
}

run_2x4x4_regs_regex <- function(panel, regex_pattern, strip_prefix, rhs, lags = NW_LAGS) {
  asset_cols <- grep(regex_pattern, names(panel), value = TRUE)
  stopifnot(length(asset_cols) == 32)

  X <- panel[, ..rhs]
  res <- list()

  for (col in asset_cols) {
    lab <- parse_2x4x4_from_col(col, strip_prefix)
    if (is.null(lab)) next

    y <- panel[[col]] - panel[["RF"]]
    est <- run_ts_reg_nw(y, X, lags)
    if (is.null(est)) next

    est[, portfolio := col]
    est[, size := lab$size]
    est[, i := lab$i]
    est[, j := lab$j]
    res[[col]] <- est
  }

  rbindlist(res, fill = TRUE)
}

run_block_model <- function(dt, asset_cols, factor_cols) {
  stopifnot("RF" %in% names(dt))
  stopifnot(all(asset_cols %in% names(dt)))
  stopifnot(all(factor_cols %in% names(dt)))

  tmp <- copy(dt[, c("RF", asset_cols, factor_cols), with = FALSE])

  for (col in asset_cols) {
    tmp[, (col) := get(col) - RF]
  }

  tmp[, RF := NULL]
  tmp <- tmp[complete.cases(tmp)]

  Tn <- nrow(tmp)
  N <- length(asset_cols)
  K <- length(factor_cols)

  if (Tn <= (N + K)) {
    return(NULL)
  }

  R <- as.matrix(tmp[, ..asset_cols])
  F <- as.matrix(tmp[, ..factor_cols])
  X <- cbind(1, F)

  XtX_inv <- solve(crossprod(X))
  B <- XtX_inv %*% crossprod(X, R)
  alpha <- as.numeric(B[1, ])

  E <- R - X %*% B
  Sigma_e <- crossprod(E) / (Tn - K - 1)
  mu_f <- colMeans(F)
  Sigma_f <- cov(F)

  if (qr(Sigma_e)$rank < ncol(Sigma_e)) {
    return(NULL)
  }
  if (qr(Sigma_f)$rank < ncol(Sigma_f)) {
    return(NULL)
  }

  alpha_mat <- matrix(alpha, ncol = 1)
  mu_f_mat <- matrix(mu_f, ncol = 1)

  grs_num <- ((Tn - N - K) / N) * t(alpha_mat) %*% solve(Sigma_e) %*% alpha_mat
  grs_den <- 1 + t(mu_f_mat) %*% solve(Sigma_f) %*% mu_f_mat

  GRS <- as.numeric(grs_num / grs_den)
  GRS_pval <- 1 - pf(GRS, N, Tn - N - K)

  abs_alpha <- abs(alpha)
  rms_alpha <- sqrt(mean(alpha^2))

  mean_excess <- colMeans(R)
  rbar <- mean_excess - mean(mean_excess)

  mean_abs_alpha <- mean(abs_alpha)
  mean_abs_alpha_over_mean_abs_rbar <- mean(abs(alpha)) / mean(abs(rbar))
  mean_alpha2_over_mean_rbar2 <- mean(alpha^2) / mean(rbar^2)

  per_asset_regs <- lapply(asset_cols, function(col) {
    reg <- run_ts_reg_nw(
      y = dt[[col]] - dt[["RF"]],
      X = dt[, ..factor_cols],
      lags = NW_LAGS
    )
    if (is.null(reg)) {
      return(data.table(
        portfolio = col,
        alpha = NA_real_,
        alpha_tstat = NA_real_,
        alpha_pval = NA_real_,
        alpha_stars = "",
        adj_r2 = NA_real_
      ))
    }

    a <- reg[term == "(Intercept)"][1]
    data.table(
      portfolio = col,
      alpha = a$estimate,
      alpha_tstat = a$tstat,
      alpha_pval = a$pval,
      alpha_stars = a$stars,
      adj_r2 = a$adj_r2
    )
  })

  per_asset_dt <- rbindlist(per_asset_regs, fill = TRUE)

  data.table(
    T = Tn,
    N = N,
    K = K,
    GRS = GRS,
    GRS_pval = GRS_pval,
    mean_abs_alpha = mean_abs_alpha,
    rms_alpha = rms_alpha,
    avg_adj_r2 = mean(per_asset_dt$adj_r2, na.rm = TRUE),
    A_abs_a_over_A_abs_rbar = mean_abs_alpha_over_mean_abs_rbar,
    A_a2_over_A_rbar2 = mean_alpha2_over_mean_rbar2
  )
}

run_block_model_with_alphas <- function(dt, asset_cols, factor_cols, lags = NW_LAGS) {
  grs <- run_block_model(dt, asset_cols, factor_cols)

  per_asset_regs <- lapply(asset_cols, function(col) {
    reg <- run_ts_reg_nw(
      y = dt[[col]] - dt[["RF"]],
      X = dt[, ..factor_cols],
      lags = lags
    )

    if (is.null(reg)) {
      return(data.table(
        portfolio = col,
        alpha = NA_real_,
        alpha_tstat = NA_real_,
        alpha_pval = NA_real_,
        alpha_stars = "",
        adj_r2 = NA_real_
      ))
    }

    a <- reg[term == "(Intercept)"][1]
    data.table(
      portfolio = col,
      alpha = a$estimate,
      alpha_tstat = a$tstat,
      alpha_pval = a$pval,
      alpha_stars = a$stars,
      adj_r2 = a$adj_r2
    )
  })

  list(
    summary = grs,
    alphas = rbindlist(per_asset_regs, fill = TRUE)
  )
}

### 10.4 Alternative Factor Sets ##############################################
smq_variants <- c(
  "SMQ_firm_VW",
  "SMQ_sic_VW",
  "SMQ_sic2_VW",
  "SMQ_firm_K",
  "SMQ_sic_K",
  "SMQ_sic2_K",
  "SMQ_firm_exp_VW",
  "SMQ_sic_exp_VW",
  "SMQ_sic2_exp_VW",
  "SMQ_firm_exp_K",
  "SMQ_sic_exp_K",
  "SMQ_sic2_exp_K"
)

alt_factor_sets <- list(
  FIRM_VW     = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_VW"),
  SIC_VW      = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_VW"),
  SIC2_VW     = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_VW"),
  FIRM_K      = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_K"),
  SIC_K       = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_K"),
  SIC2_K      = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_K"),
  FIRM_EXP_VW = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_VW"),
  SIC_EXP_VW  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_VW"),
  SIC2_EXP_VW = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_VW"),
  FIRM_EXP_K  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_K"),
  SIC_EXP_K   = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_K"),
  SIC2_EXP_K  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_K")
)

alt_summary_tables <- list()
alt_corr_tables <- list()
alt_factor_on_factor <- list()

for (nm in names(alt_factor_sets)) {
  facs <- alt_factor_sets[[nm]]
  alt_summary_tables[[nm]] <- make_summary_table(panel, facs)
  alt_corr_tables[[nm]] <- cor(panel[, ..facs], use = "pairwise.complete.obs")
  alt_factor_on_factor[[nm]] <- make_factor_table(panel, facs, NW_LAGS)
}

saveRDS(alt_summary_tables, "section10_alt_factor_summary_tables.rds")
saveRDS(alt_corr_tables, "section10_alt_factor_correlations.rds")
saveRDS(alt_factor_on_factor, "section10_alt_factor_on_factor_tables.rds")

### 10.5 SMQ Comparison Across Specifications #################################
smq_summary_table <- make_smq_comparison_table(panel, smq_variants, NW_LAGS)
smq_summary_display <- make_smq_display_table(smq_summary_table)

saveRDS(smq_summary_table, "section10_smq_summary_table.rds")
saveRDS(smq_summary_display, "section10_smq_summary_display_table.rds")

### 10.6 Alternative Factor ACF ###############################################
plot_in_pages <- function(var_names, nrow = 3, ncol = 2, mar = c(3, 4, 2, 1), oma = c(0, 0, 2, 0)) {
  n_per_page <- nrow * ncol
  n_vars <- length(var_names)
  n_pages <- ceiling(n_vars / n_per_page)

  for (pg in seq_len(n_pages)) {
    idx_start <- (pg - 1) * n_per_page + 1
    idx_end <- min(pg * n_per_page, n_vars)
    vars_pg <- var_names[idx_start:idx_end]

    par(mfrow = c(nrow, ncol), mar = mar, oma = oma)

    for (f in vars_pg) {
      acf(panel[[f]], na.action = na.pass, main = paste("ACF:", f))
    }

    n_empty <- n_per_page - length(vars_pg)
    if (n_empty > 0) {
      for (i in seq_len(n_empty)) {
        plot.new()
      }
    }

    mtext(
      text = paste("Alternative Factor ACFs - Page", pg, "of", n_pages),
      outer = TRUE,
      cex = 1
    )
  }

  par(mfrow = c(1, 1))
}

alt_factor_cols <- smq_variants[smq_variants %in% names(panel)]
plot_in_pages(alt_factor_cols, nrow = 3, ncol = 2)

### 10.7 Alternative Factor Time-Series Plots #################################
par(mfrow = c(ceiling(length(alt_factor_cols) / 2), 2), mar = c(3, 4, 2, 1))
for (f in alt_factor_cols) {
  plot(
    panel$mdate, panel[[f]],
    type = "l", lwd = 1.3,
    main = f, ylab = "Return", xlab = ""
  )
}
par(mfrow = c(1, 1))

### 10.8 Alternative Factor Cumulative Return Plots ###########################
par(mfrow = c(ceiling(length(alt_factor_cols) / 2), 2), mar = c(3, 4, 2, 1))
for (f in alt_factor_cols) {
  plot(
    panel$mdate, cum_safe(panel[[f]]),
    type = "l", lwd = 1.3,
    main = paste("Cumulative", f), ylab = "Cumulative Return", xlab = ""
  )
}
par(mfrow = c(1, 1))

### 10.9 Alternative Decile Strategy Tables ###################################
models_basic <- list(
  CAPM = c("MKT"),
  FF3  = c("MKT", "SMB", "HML"),
  FF5  = c("MKT", "SMB", "HML", "RMW", "CMA")
)

build_decile_table_with_smq <- function(panel, dec_prefix, smq_col, lags = NW_LAGS, digits_est = 2, digits_t = 2) {
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
    vol_cells[j] <- sprintf("%.2f", 100 * sd_m * sqrt(12))
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
  vol_cells[11] <- sprintf("%.2f", 100 * sd_smq * sqrt(12))
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

  capm_cells <- alpha_row(models_basic$CAPM)
  ff3_cells <- alpha_row(models_basic$FF3)
  ff5_cells <- alpha_row(models_basic$FF5)

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

alt_decile_tables <- list(
  CHI_FIRM_DEC10_VW     = build_decile_table_with_smq(factor_test_panel, "CHI_FIRM_DD", "SMQ_firm_VW", NW_LAGS),
  CHI_SIC_DEC10_VW      = build_decile_table_with_smq(factor_test_panel, "CHI_SIC_DD", "SMQ_sic_VW", NW_LAGS),
  CHI_SIC2_DEC10_VW     = build_decile_table_with_smq(factor_test_panel, "CHI_SIC2_DD", "SMQ_sic2_VW", NW_LAGS),
  CHI_FIRM_EXP_DEC10_VW = build_decile_table_with_smq(factor_test_panel, "CHI_FIRM_EXP_DD", "SMQ_firm_exp_VW", NW_LAGS),
  CHI_SIC_EXP_DEC10_VW  = build_decile_table_with_smq(factor_test_panel, "CHI_SIC_EXP_DD", "SMQ_sic_exp_VW", NW_LAGS),
  CHI_SIC2_EXP_DEC10_VW = build_decile_table_with_smq(factor_test_panel, "CHI_SIC2_EXP_DD", "SMQ_sic2_exp_VW", NW_LAGS)
)

saveRDS(alt_decile_tables, "section10_alt_decile_tables.rds")

### 10.10 Alternative Decile Cell Counts ######################################
decile_vars <- c(
  "chi_firm_dec10",
  "chi_sic_dec10",
  "chi_sic2_dec10",
  "chi_firm_exp_dec10",
  "chi_sic_exp_dec10",
  "chi_sic2_exp_dec10"
)

count_tables_deciles <- list()

for (v in decile_vars) {
  tmp <- labelled_data[
    !is.na(get(v)),
    .(N = uniqueN(permno)),
    by = .(mdate, decile = get(v))
  ][
    ,
    .(Avg_N = mean(N, na.rm = TRUE)),
    by = decile
  ][order(decile)]

  count_tables_deciles[[v]] <- tmp
}

saveRDS(count_tables_deciles, "section10_alt_decile_cell_counts.rds")

### 10.11 Alternative Leg Diagnostics #########################################
june_form <- raw_data[
  month == 6 &
    has_dec == TRUE &
    has_june == TRUE,
  .(
    permno,
    gvkey,
    sic,
    sic2,
    ffyear        = ffyear + 1L,
    exchcd,
    me_june       = me_clean,
    bm            = bm_ff,
    op            = op_clean,
    inv           = inv_annual,
    ik            = ik_annual,
    q             = tobins_q,
    chi_firm      = chi_firm,
    chi_firm_exp  = chi_firm_exp,
    chi_sic       = chi_sic,
    chi_sic_exp   = chi_sic_exp,
    chi_sic2      = chi_sic2,
    chi_sic2_exp  = chi_sic2_exp
  )
]

signals_nyse <- june_form[exchcd == 1]

chi_firm_bp <- signals_nyse[
  !is.na(chi_firm),
  {
    b <- ff_30_70_breaks(chi_firm)
    .(chi_firm_30 = b["low"], chi_firm_70 = b["high"])
  },
  by = ffyear
]

chi_firm_exp_bp <- signals_nyse[
  !is.na(chi_firm_exp),
  {
    b <- ff_30_70_breaks(chi_firm_exp)
    .(chi_firm_exp_30 = b["low"], chi_firm_exp_70 = b["high"])
  },
  by = ffyear
]

chi_sic_bp <- signals_nyse[
  !is.na(chi_sic),
  {
    b <- ff_30_70_breaks(chi_sic)
    .(chi_sic_30 = b["low"], chi_sic_70 = b["high"])
  },
  by = ffyear
]

chi_sic_exp_bp <- signals_nyse[
  !is.na(chi_sic_exp),
  {
    b <- ff_30_70_breaks(chi_sic_exp)
    .(chi_sic_exp_30 = b["low"], chi_sic_exp_70 = b["high"])
  },
  by = ffyear
]

chi_sic2_bp <- signals_nyse[
  !is.na(chi_sic2),
  {
    b <- ff_30_70_breaks(chi_sic2)
    .(chi_sic2_30 = b["low"], chi_sic2_70 = b["high"])
  },
  by = ffyear
]

chi_sic2_exp_bp <- signals_nyse[
  !is.na(chi_sic2_exp),
  {
    b <- ff_30_70_breaks(chi_sic2_exp)
    .(chi_sic2_exp_30 = b["low"], chi_sic2_exp_70 = b["high"])
  },
  by = ffyear
]

june_form <- Reduce(
  function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
  list(
    june_form,
    chi_firm_bp,
    chi_firm_exp_bp,
    chi_sic_bp,
    chi_sic_exp_bp,
    chi_sic2_bp,
    chi_sic2_exp_bp
  )
)

june_form[, chi_firm_grp := fcase(
  !is.na(chi_firm) & chi_firm <= chi_firm_30, "Low",
  !is.na(chi_firm) & chi_firm > chi_firm_70, "High",
  default = NA_character_
)]

june_form[, chi_firm_exp_grp := fcase(
  !is.na(chi_firm_exp) & chi_firm_exp <= chi_firm_exp_30, "Low",
  !is.na(chi_firm_exp) & chi_firm_exp > chi_firm_exp_70, "High",
  default = NA_character_
)]

june_form[, chi_sic_grp := fcase(
  !is.na(chi_sic) & chi_sic <= chi_sic_30, "Low",
  !is.na(chi_sic) & chi_sic > chi_sic_70, "High",
  default = NA_character_
)]

june_form[, chi_sic_exp_grp := fcase(
  !is.na(chi_sic_exp) & chi_sic_exp <= chi_sic_exp_30, "Low",
  !is.na(chi_sic_exp) & chi_sic_exp > chi_sic_exp_70, "High",
  default = NA_character_
)]

june_form[, chi_sic2_grp := fcase(
  !is.na(chi_sic2) & chi_sic2 <= chi_sic2_30, "Low",
  !is.na(chi_sic2) & chi_sic2 > chi_sic2_70, "High",
  default = NA_character_
)]

june_form[, chi_sic2_exp_grp := fcase(
  !is.na(chi_sic2_exp) & chi_sic2_exp <= chi_sic2_exp_30, "Low",
  !is.na(chi_sic2_exp) & chi_sic2_exp > chi_sic2_exp_70, "High",
  default = NA_character_
)]

### 10.11A Alternative Leg Industry Concentration: Count-Based HHI ############
factor_leg_hhi_alt <- rbindlist(list(
  make_leg_hhi_table(june_form, "chi_firm_grp", "SMQ_firm", "sic2"),
  make_leg_hhi_table(june_form, "chi_firm_exp_grp", "SMQ_firm_exp", "sic2"),
  make_leg_hhi_table(june_form, "chi_sic_grp", "SMQ_sic", "sic2"),
  make_leg_hhi_table(june_form, "chi_sic_exp_grp", "SMQ_sic_exp", "sic2"),
  make_leg_hhi_table(june_form, "chi_sic2_grp", "SMQ_sic2", "sic2"),
  make_leg_hhi_table(june_form, "chi_sic2_exp_grp", "SMQ_sic2_exp", "sic2")
), fill = TRUE)

factor_leg_hhi_summary_alt <- make_leg_hhi_summary(
  factor_leg_hhi_alt,
  hhi_col = "HHI",
  effn_col = "Effective_N"
)

saveRDS(factor_leg_hhi_alt, "section10_alt_factor_leg_hhi_by_year.rds")
saveRDS(factor_leg_hhi_summary_alt, "section10_alt_factor_leg_hhi_summary.rds")

### 10.11B Alternative Leg Industry Concentration: Value-Weighted HHI #########
factor_leg_hhi_weighted_alt <- rbindlist(list(
  make_leg_weighted_hhi_table(june_form, "chi_firm_grp", "SMQ_firm", "sic2", "me_june"),
  make_leg_weighted_hhi_table(june_form, "chi_firm_exp_grp", "SMQ_firm_exp", "sic2", "me_june"),
  make_leg_weighted_hhi_table(june_form, "chi_sic_grp", "SMQ_sic", "sic2", "me_june"),
  make_leg_weighted_hhi_table(june_form, "chi_sic_exp_grp", "SMQ_sic_exp", "sic2", "me_june"),
  make_leg_weighted_hhi_table(june_form, "chi_sic2_grp", "SMQ_sic2", "sic2", "me_june"),
  make_leg_weighted_hhi_table(june_form, "chi_sic2_exp_grp", "SMQ_sic2_exp", "sic2", "me_june")
), fill = TRUE)

factor_leg_hhi_weighted_summary_alt <- make_leg_hhi_summary(
  factor_leg_hhi_weighted_alt,
  hhi_col = "HHI_weighted",
  effn_col = "Effective_N_weighted"
)

saveRDS(
  factor_leg_hhi_weighted_alt,
  "section10_alt_factor_leg_hhi_weighted_by_year.rds"
)
saveRDS(
  factor_leg_hhi_weighted_summary_alt,
  "section10_alt_factor_leg_hhi_weighted_summary.rds"
)

### 10.11C Dominant industries within each leg ################################
make_leg_industry_shares <- function(dt, grp_var, factor_name, industry_var = "sic2", weight_var = "me_june") {
  tmp <- dt[
    !is.na(get(grp_var)) &
      !is.na(get(industry_var)) &
      !is.na(get(weight_var)) &
      is.finite(get(weight_var)) &
      get(weight_var) > 0,
    .(
      industry_weight = sum(get(weight_var), na.rm = TRUE)
    ),
    by = .(ffyear, Leg = get(grp_var), industry = get(industry_var))
  ]

  tmp[, total_leg_weight := sum(industry_weight), by = .(ffyear, Leg)]
  tmp[, share := industry_weight / total_leg_weight]
  tmp[, hhi_component := share^2]
  tmp[, rank_in_leg := frank(-share, ties.method = "min"), by = .(ffyear, Leg)]
  tmp[, Factor := factor_name]

  setcolorder(
    tmp,
    c("Factor", "ffyear", "Leg", "industry", "industry_weight", "total_leg_weight", "share", "hhi_component", "rank_in_leg")
  )

  tmp[]
}

make_leg_industry_summary <- function(industry_share_dt) {
  industry_share_dt[
    ,
    .(
      industry_name = {
        vals <- unique(na.omit(industry_name))
        if (length(vals) == 0) NA_character_ else vals[1]
      },
      Avg_share = mean(share, na.rm = TRUE),
      Median_share = median(share, na.rm = TRUE),
      Max_share = max(share, na.rm = TRUE),
      Avg_HHI_component = mean(hhi_component, na.rm = TRUE),
      Times_top1 = sum(rank_in_leg == 1, na.rm = TRUE),
      Share_top1_years = mean(rank_in_leg == 1, na.rm = TRUE),
      Avg_rank = mean(rank_in_leg, na.rm = TRUE),
      N_years_present = uniqueN(ffyear)
    ),
    by = .(Factor, Leg, industry)
  ][order(Factor, Leg, -Avg_share, -Times_top1)]
}

get_top_n_industries <- function(industry_summary_dt, n = 10) {
  industry_summary_dt[
    order(Factor, Leg, -Avg_share, -Times_top1),
    head(.SD, n),
    by = .(Factor, Leg)
  ]
}

make_topk_share_summary <- function(industry_share_dt, k = 3) {
  industry_share_dt[
    rank_in_leg <= k,
    .(
      TopK_share = sum(share, na.rm = TRUE)
    ),
    by = .(Factor, ffyear, Leg)
  ][
    ,
    .(
      Avg_TopK_share = mean(TopK_share, na.rm = TRUE),
      Median_TopK_share = median(TopK_share, na.rm = TRUE),
      Max_TopK_share = max(TopK_share, na.rm = TRUE)
    ),
    by = .(Factor, Leg)
  ][order(Factor, Leg)]
}

sic2_lookup <- data.table(
  industry = c(
    1, 2, 7, 8, 9,
    10, 12, 13, 14,
    20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
    30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
    40, 41, 42, 44, 45, 46, 47, 48, 49,
    50, 51,
    52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 67,
    70, 72, 73, 75, 76, 78, 79,
    80, 81, 82, 83, 84, 86, 87, 89,
    91, 92, 93, 94, 95, 96, 97, 99
  ),
  industry_name = c(
    "Agricultural production - crops",
    "Agricultural production - livestock",
    "Agricultural services",
    "Forestry",
    "Fishing, hunting, trapping",
    "Metal mining",
    "Coal mining",
    "Oil and gas extraction",
    "Nonmetallic minerals",
    "Food products",
    "Tobacco products",
    "Textile mill products",
    "Apparel",
    "Lumber and wood",
    "Furniture",
    "Paper",
    "Printing and publishing",
    "Chemicals",
    "Petroleum refining",
    "Rubber and plastics",
    "Leather",
    "Stone, clay, glass",
    "Primary metals",
    "Fabricated metal products",
    "Industrial machinery",
    "Electronic equipment",
    "Transportation equipment",
    "Instruments",
    "Miscellaneous manufacturing",
    "Railroad transportation",
    "Local passenger transit",
    "Motor freight",
    "Water transportation",
    "Air transportation",
    "Pipelines",
    "Transportation services",
    "Communications",
    "Utilities",
    "Wholesale - durable goods",
    "Wholesale - nondurable goods",
    "Retail - building materials",
    "Retail - general merchandise",
    "Retail - food stores",
    "Retail - automotive dealers",
    "Retail - apparel",
    "Retail - home furniture",
    "Retail - eating and drinking",
    "Retail - miscellaneous",
    "Depository institutions",
    "Nondepository institutions",
    "Security and commodity brokers",
    "Insurance carriers",
    "Insurance agents",
    "Real estate",
    "Holding and other investment offices",
    "Hotels and lodging",
    "Personal services",
    "Business services",
    "Auto repair and parking",
    "Miscellaneous repair",
    "Motion pictures",
    "Amusement and recreation",
    "Health services",
    "Legal services",
    "Educational services",
    "Social services",
    "Museums and membership organizations",
    "Engineering and management services",
    "Private households",
    "Executive, legislative, general government",
    "Justice, public order, safety",
    "Public finance, taxation",
    "Administration of human resources",
    "Environmental quality and housing",
    "Administration of economic programs",
    "National security and international affairs",
    "Nonclassifiable establishments"
  )
)

factor_leg_industry_shares_alt <- rbindlist(list(
  make_leg_industry_shares(june_form, "chi_firm_grp", "SMQ_firm", "sic2", "me_june"),
  make_leg_industry_shares(june_form, "chi_firm_exp_grp", "SMQ_firm_exp", "sic2", "me_june"),
  make_leg_industry_shares(june_form, "chi_sic_grp", "SMQ_sic", "sic2", "me_june"),
  make_leg_industry_shares(june_form, "chi_sic_exp_grp", "SMQ_sic_exp", "sic2", "me_june"),
  make_leg_industry_shares(june_form, "chi_sic2_grp", "SMQ_sic2", "sic2", "me_june"),
  make_leg_industry_shares(june_form, "chi_sic2_exp_grp", "SMQ_sic2_exp", "sic2", "me_june")
), fill = TRUE)

factor_leg_industry_shares_alt <- merge(
  factor_leg_industry_shares_alt,
  sic2_lookup,
  by = "industry",
  all.x = TRUE
)

factor_leg_industry_summary_alt <- make_leg_industry_summary(factor_leg_industry_shares_alt)

factor_leg_top10_industries_alt <- get_top_n_industries(
  factor_leg_industry_summary_alt,
  n = 10
)

factor_leg_top1_share_alt <- make_topk_share_summary(
  factor_leg_industry_shares_alt,
  k = 1
)

factor_leg_top3_share_alt <- make_topk_share_summary(
  factor_leg_industry_shares_alt,
  k = 3
)

factor_leg_dominant_industry_by_year_alt <- factor_leg_industry_shares_alt[
  rank_in_leg == 1
][order(Factor, ffyear, Leg)]

saveRDS(
  factor_leg_industry_shares_alt,
  "section10_alt_factor_leg_industry_shares_by_year.rds"
)

saveRDS(
  factor_leg_industry_summary_alt,
  "section10_alt_factor_leg_industry_summary.rds"
)

saveRDS(
  factor_leg_top10_industries_alt,
  "section10_alt_factor_leg_top10_industries.rds"
)

saveRDS(
  factor_leg_top1_share_alt,
  "section10_alt_factor_leg_top1_share_summary.rds"
)

saveRDS(
  factor_leg_top3_share_alt,
  "section10_alt_factor_leg_top3_share_summary.rds"
)

saveRDS(
  factor_leg_dominant_industry_by_year_alt,
  "section10_alt_factor_leg_dominant_industry_by_year.rds"
)

factor_leg_counts_alt <- rbindlist(list(
  june_form[!is.na(chi_firm_grp), .(N = .N), by = .(ffyear, Leg = chi_firm_grp)][
    , .(Factor = "SMQ_firm", Leg, Avg_N = mean(N), Median_N = median(N))
  ],
  june_form[!is.na(chi_firm_exp_grp), .(N = .N), by = .(ffyear, Leg = chi_firm_exp_grp)][
    , .(Factor = "SMQ_firm_exp", Leg, Avg_N = mean(N), Median_N = median(N))
  ],
  june_form[!is.na(chi_sic_grp), .(N = .N), by = .(ffyear, Leg = chi_sic_grp)][
    , .(Factor = "SMQ_sic", Leg, Avg_N = mean(N), Median_N = median(N))
  ],
  june_form[!is.na(chi_sic_exp_grp), .(N = .N), by = .(ffyear, Leg = chi_sic_exp_grp)][
    , .(Factor = "SMQ_sic_exp", Leg, Avg_N = mean(N), Median_N = median(N))
  ],
  june_form[!is.na(chi_sic2_grp), .(N = .N), by = .(ffyear, Leg = chi_sic2_grp)][
    , .(Factor = "SMQ_sic2", Leg, Avg_N = mean(N), Median_N = median(N))
  ],
  june_form[!is.na(chi_sic2_exp_grp), .(N = .N), by = .(ffyear, Leg = chi_sic2_exp_grp)][
    , .(Factor = "SMQ_sic2_exp", Leg, Avg_N = mean(N), Median_N = median(N))
  ]
), fill = TRUE)

saveRDS(factor_leg_counts_alt, "section10_alt_factor_leg_counts.rds")

leg_spread_table_alt <- rbindlist(list(
  make_leg_spread_table(june_form, "chi_firm_grp", "SMQ_firm", "chi_firm"),
  make_leg_spread_table(june_form, "chi_firm_exp_grp", "SMQ_firm_exp", "chi_firm_exp"),
  make_leg_spread_table(june_form, "chi_sic_grp", "SMQ_sic", "chi_sic"),
  make_leg_spread_table(june_form, "chi_sic_exp_grp", "SMQ_sic_exp", "chi_sic_exp"),
  make_leg_spread_table(june_form, "chi_sic2_grp", "SMQ_sic2", "chi_sic2"),
  make_leg_spread_table(june_form, "chi_sic2_exp_grp", "SMQ_sic2_exp", "chi_sic2_exp")
), fill = TRUE)

saveRDS(leg_spread_table_alt, "section10_alt_leg_spread_table.rds")

### 10.12 Alternative Leg Return Statistics ###################################
leg_labels <- june_form[, .(
  permno,
  ffyear,
  chi_firm_grp,
  chi_firm_exp_grp,
  chi_sic_grp,
  chi_sic_exp_grp,
  chi_sic2_grp,
  chi_sic2_exp_grp
)]

monthly_alt <- merge(
  raw_data,
  leg_labels,
  by = c("permno", "ffyear"),
  all.x = TRUE
)

setorder(monthly_alt, permno, mdate)
monthly_alt[, me_ff_w := shift(me_clean), by = permno]
monthly_alt[me_ff_w <= 0, me_ff_w := NA_real_]

tmp_firm <- make_leg_return_stats(monthly_alt, "chi_firm_grp", "SMQ_firm")
tmp_firm_exp <- make_leg_return_stats(monthly_alt, "chi_firm_exp_grp", "SMQ_firm_exp")
tmp_sic <- make_leg_return_stats(monthly_alt, "chi_sic_grp", "SMQ_sic")
tmp_sic_exp <- make_leg_return_stats(monthly_alt, "chi_sic_exp_grp", "SMQ_sic_exp")
tmp_sic2 <- make_leg_return_stats(monthly_alt, "chi_sic2_grp", "SMQ_sic2")
tmp_sic2_exp <- make_leg_return_stats(monthly_alt, "chi_sic2_exp_grp", "SMQ_sic2_exp")

factor_leg_return_stats_alt <- rbindlist(list(
  tmp_firm$stats,
  tmp_firm_exp$stats,
  tmp_sic$stats,
  tmp_sic_exp$stats,
  tmp_sic2$stats,
  tmp_sic2_exp$stats
), fill = TRUE)

spread_series_alt <- list(
  SMQ_firm     = tmp_firm$spread_series,
  SMQ_firm_exp = tmp_firm_exp$spread_series,
  SMQ_sic      = tmp_sic$spread_series,
  SMQ_sic_exp  = tmp_sic_exp$spread_series,
  SMQ_sic2     = tmp_sic2$spread_series,
  SMQ_sic2_exp = tmp_sic2_exp$spread_series
)

saveRDS(factor_leg_return_stats_alt, "section10_alt_factor_leg_return_stats.rds")
saveRDS(spread_series_alt, "section10_alt_spread_series.rds")

### 10.13 Alternative Spread Histograms and Densities #########################
plot_names_hist <- c("SMQ_firm", "SMQ_firm_exp", "SMQ_sic", "SMQ_sic_exp", "SMQ_sic2", "SMQ_sic2_exp")

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (nm in plot_names_hist) {
  hist(
    100 * spread_series_alt[[nm]]$exret,
    breaks = 30,
    main = paste0(nm, "\nHigh-Low"),
    xlab = "Monthly excess return (%)"
  )
}
par(mfrow = c(1, 1))

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
for (nm in plot_names_hist) {
  x <- 100 * spread_series_alt[[nm]]$exret
  plot(density(x, na.rm = TRUE), main = paste("Density:", nm), xlab = "Monthly excess return (%)")
  abline(v = mean(x, na.rm = TRUE), lty = 2)
}
par(mfrow = c(1, 1))

### 10.14 Alternative 5x5 Portfolio Performance ###############################
port_5x5_specs_alt <- list(
  chi_firm = list(regex = "^CHI_FIRM_S[1-5]C[1-5]$", prefix = "CHI_FIRM_"),
  chi_firm_exp = list(regex = "^CHI_FIRM_EXP_S[1-5]C[1-5]$", prefix = "CHI_FIRM_EXP_"),
  chi_sic = list(regex = "^CHI_SIC_S[1-5]C[1-5]$", prefix = "CHI_SIC_"),
  chi_sic_exp = list(regex = "^CHI_SIC_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC_EXP_"),
  chi_sic2 = list(regex = "^CHI_SIC2_S[1-5]C[1-5]$", prefix = "CHI_SIC2_"),
  chi_sic2_exp = list(regex = "^CHI_SIC2_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC2_EXP_")
)

port_5x5_alt <- lapply(port_5x5_specs_alt, function(sp) {
  wide_5x5_to_long_regex(factor_test_panel, sp$regex, sp$prefix)
})

mean_5x5_alt <- list()
mean_excess_5x5_alt <- list()
sd_5x5_alt <- list()
sr_5x5_alt <- list()

for (nm in names(port_5x5_alt)) {
  mean_5x5_alt[[nm]] <- make_5x5_mean_matrix(filter_window(port_5x5_alt[[nm]]), "size5", "char5")
  tmp <- build_excess_stats(port_5x5_alt[[nm]], "char5", rf_monthly)
  mean_excess_5x5_alt[[nm]] <- tmp$mean
  sd_5x5_alt[[nm]] <- tmp$sd
  sr_5x5_alt[[nm]] <- tmp$sr
}

saveRDS(mean_5x5_alt, "section10_mean_5x5_alt.rds")
saveRDS(mean_excess_5x5_alt, "section10_mean_excess_5x5_alt.rds")
saveRDS(sd_5x5_alt, "section10_sd_5x5_alt.rds")
saveRDS(sr_5x5_alt, "section10_sr_5x5_alt.rds")

### 10.15 HMLO Construction ###################################################
align_orthogonalized_return <- function(fit, n_rows) {
  values <- rep(NA_real_, n_rows)
  used_rows <- as.integer(names(residuals(fit)))
  values[used_rows] <- unname(coef(fit)[1]) + residuals(fit)
  values
}

hmlo_fit_firm <- lm(HML ~ MKT + SMB + RMW + CMA + SMQ_firm_VW, data = factor_test_panel)
factor_test_panel[, HMLO_FIRM := align_orthogonalized_return(hmlo_fit_firm, .N)]

hmlo_fit_sic <- lm(HML ~ MKT + SMB + RMW + CMA + SMQ_sic_VW, data = factor_test_panel)
factor_test_panel[, HMLO_SIC := align_orthogonalized_return(hmlo_fit_sic, .N)]

hmlo_fit_sic2 <- lm(HML ~ MKT + SMB + RMW + CMA + SMQ_sic2_VW, data = factor_test_panel)
factor_test_panel[, HMLO_SIC2 := align_orthogonalized_return(hmlo_fit_sic2, .N)]

### 10.16 Benchmark Model Definitions #########################################
models_general <- list(
  CAPM = c("MKT"),
  FF3  = c("MKT", "SMB", "HML"),
  FF5  = c("MKT", "SMB", "HML", "RMW", "CMA")
)

models_5x5_alt <- list(
  FF5_SMQ_FIRM_VW     = c("MKT", "SMB", "HMLO_FIRM", "RMW", "CMA", "SMQ_firm_VW"),
  FF5_SMQ_SIC_VW      = c("MKT", "SMB", "HMLO_SIC", "RMW", "CMA", "SMQ_sic_VW"),
  FF5_SMQ_SIC2_VW     = c("MKT", "SMB", "HMLO_SIC2", "RMW", "CMA", "SMQ_sic2_VW"),
  FF5_SMQ_FIRM_K      = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_K"),
  FF5_SMQ_SIC_K       = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_K"),
  FF5_SMQ_SIC2_K      = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_K"),
  FF5_SMQ_FIRM_EXP_VW = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_VW"),
  FF5_SMQ_SIC_EXP_VW  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_VW"),
  FF5_SMQ_SIC2_EXP_VW = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_VW"),
  FF5_SMQ_FIRM_EXP_K  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_K"),
  FF5_SMQ_SIC_EXP_K   = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_K"),
  FF5_SMQ_SIC2_EXP_K  = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_K")
)

variant_factor_map <- list(
  FIRM_VW     = list(smq_col = "SMQ_firm_VW", ff5_smq = c("MKT", "SMB", "HMLO_FIRM", "RMW", "CMA", "SMQ_firm_VW")),
  SIC_VW      = list(smq_col = "SMQ_sic_VW", ff5_smq = c("MKT", "SMB", "HMLO_SIC", "RMW", "CMA", "SMQ_sic_VW")),
  SIC2_VW     = list(smq_col = "SMQ_sic2_VW", ff5_smq = c("MKT", "SMB", "HMLO_SIC2", "RMW", "CMA", "SMQ_sic2_VW")),
  FIRM_K      = list(smq_col = "SMQ_firm_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_K")),
  SIC_K       = list(smq_col = "SMQ_sic_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_K")),
  SIC2_K      = list(smq_col = "SMQ_sic2_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_K")),
  FIRM_EXP_VW = list(smq_col = "SMQ_firm_exp_VW", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_VW")),
  SIC_EXP_VW  = list(smq_col = "SMQ_sic_exp_VW", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_VW")),
  SIC2_EXP_VW = list(smq_col = "SMQ_sic2_exp_VW", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_VW")),
  FIRM_EXP_K  = list(smq_col = "SMQ_firm_exp_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_firm_exp_K")),
  SIC_EXP_K   = list(smq_col = "SMQ_sic_exp_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic_exp_K")),
  SIC2_EXP_K  = list(smq_col = "SMQ_sic2_exp_K", ff5_smq = c("MKT", "SMB", "HML", "RMW", "CMA", "SMQ_sic2_exp_K"))
)

### 10.17 Alternative 5x5 Alpha Tables ########################################
asset_blocks_5x5_alt <- list(
  BM            = list(regex = "^BM_S[1-5]V[1-5]$", prefix = "BM_"),
  OP            = list(regex = "^OP_S[1-5]P[1-5]$", prefix = "OP_"),
  INV           = list(regex = "^INV_S[1-5]I[1-5]$", prefix = "INV_"),
  CHI_FIRM      = list(regex = "^CHI_FIRM_S[1-5]C[1-5]$", prefix = "CHI_FIRM_"),
  CHI_FIRM_EXP  = list(regex = "^CHI_FIRM_EXP_S[1-5]C[1-5]$", prefix = "CHI_FIRM_EXP_"),
  CHI_SIC       = list(regex = "^CHI_SIC_S[1-5]C[1-5]$", prefix = "CHI_SIC_"),
  CHI_SIC_EXP   = list(regex = "^CHI_SIC_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC_EXP_"),
  CHI_SIC2      = list(regex = "^CHI_SIC2_S[1-5]C[1-5]$", prefix = "CHI_SIC2_"),
  CHI_SIC2_EXP  = list(regex = "^CHI_SIC2_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC2_EXP_")
)

pick_5x5_robust_models <- function(asset_class) {
  c(
    "FF5_SMQ_FIRM_VW",
    "FF5_SMQ_SIC_VW",
    "FF5_SMQ_SIC2_VW",
    "FF5_SMQ_FIRM_K",
    "FF5_SMQ_SIC_K",
    "FF5_SMQ_SIC2_K",
    "FF5_SMQ_FIRM_EXP_VW",
    "FF5_SMQ_SIC_EXP_VW",
    "FF5_SMQ_SIC2_EXP_VW",
    "FF5_SMQ_FIRM_EXP_K",
    "FF5_SMQ_SIC_EXP_K",
    "FF5_SMQ_SIC2_EXP_K"
  )
}

robustness_outputs_5x5 <- list()
robustness_alpha_regs_5x5 <- list()

for (fam in names(asset_blocks_5x5_alt)) {
  sp <- asset_blocks_5x5_alt[[fam]]
  out <- list()
  out_reg <- list()

  for (model_name in pick_5x5_robust_models(fam)) {
    reg <- run_asset_block_regs_regex(
      factor_test_panel,
      sp$regex,
      sp$prefix,
      models_5x5_alt[[model_name]],
      NW_LAGS
    )

    out[[model_name]] <- term_to_5x5_matrix(
      reg,
      "(Intercept)",
      digits_est = 2,
      digits_t = 2,
      scale_est = 100
    )

    out_reg[[model_name]] <- extract_alpha_table(reg)
  }

  robustness_outputs_5x5[[fam]] <- out
  robustness_alpha_regs_5x5[[fam]] <- out_reg
}

saveRDS(robustness_outputs_5x5, "section10_robustness_intercepts_5x5.rds")
saveRDS(robustness_alpha_regs_5x5, "section10_robustness_intercepts_5x5_alphas.rds")

### 10.18 Print Alternative 5x5 Performance ###################################
display_all_5x5(mean_5x5_alt, "SECTION 10: MONTHLY MEAN 5x5 RETURNS (%)")
display_all_5x5(mean_excess_5x5_alt, "SECTION 10: MONTHLY MEAN 5x5 EXCESS RETURNS (%)")
display_all_5x5(sd_5x5_alt, "SECTION 10: MONTHLY STD DEV 5x5 EXCESS RETURNS (%)")
display_all_5x5(sr_5x5_alt, "SECTION 10: SHARPE RATIOS 5x5 (monthly)")

### 10.19 Alternative 2x4x4 Portfolio Performance #############################
families_2x4x4_alt_desc <- list(
  OP_CHI_FIRM  = list(regex = "^OP_CHI_FIRM_[SB][1-4][1-4]$", prefix = "OP_CHI_FIRM_"),
  OP_CHI_SIC   = list(regex = "^OP_CHI_SIC_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC_"),
  OP_CHI_SIC2  = list(regex = "^OP_CHI_SIC2_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC2_"),
  BM_CHI_FIRM  = list(regex = "^BM_CHI_FIRM_[SB][1-4][1-4]$", prefix = "BM_CHI_FIRM_"),
  BM_CHI_SIC   = list(regex = "^BM_CHI_SIC_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC_"),
  BM_CHI_SIC2  = list(regex = "^BM_CHI_SIC2_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC2_"),
  INV_CHI_FIRM = list(regex = "^INV_CHI_FIRM_[SB][1-4][1-4]$", prefix = "INV_CHI_FIRM_"),
  INV_CHI_SIC  = list(regex = "^INV_CHI_SIC_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC_"),
  INV_CHI_SIC2 = list(regex = "^INV_CHI_SIC2_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC2_")
)

long_2x4x4_alt_desc <- setNames(
  lapply(names(families_2x4x4_alt_desc), function(nm) {
    sp <- families_2x4x4_alt_desc[[nm]]
    wide_2x4x4_to_long_regex(factor_test_panel, sp$regex, sp$prefix)
  }),
  names(families_2x4x4_alt_desc)
)

long_2x4x4_alt_desc <- long_2x4x4_alt_desc[!vapply(long_2x4x4_alt_desc, is.null, logical(1))]

mean_2x4x4_alt <- list()
mean_excess_2x4x4_alt <- list()
sd_2x4x4_alt <- list()
sr_2x4x4_alt <- list()

for (nm in names(long_2x4x4_alt_desc)) {
  mean_2x4x4_alt[[nm]] <- make_mean_table_2x4x4(long_2x4x4_alt_desc[[nm]])
  mean_excess_2x4x4_alt[[nm]] <- make_mean_excess_table_2x4x4(long_2x4x4_alt_desc[[nm]], rf_monthly)
  sd_2x4x4_alt[[nm]] <- make_sd_table_2x4x4(long_2x4x4_alt_desc[[nm]], rf_monthly)
  sr_2x4x4_alt[[nm]] <- make_sr_table_2x4x4(long_2x4x4_alt_desc[[nm]], rf_monthly)
}

saveRDS(mean_2x4x4_alt, "section10_mean_2x4x4_alt.rds")
saveRDS(mean_excess_2x4x4_alt, "section10_mean_excess_2x4x4_alt.rds")
saveRDS(sd_2x4x4_alt, "section10_sd_2x4x4_alt.rds")
saveRDS(sr_2x4x4_alt, "section10_sr_2x4x4_alt.rds")

### 10.20 Alternative 2x4x4 Alpha Tables ######################################
models_2x4x4_alt <- models_5x5_alt

asset_blocks_2x4x4_alt <- list(
  BM_OP = list(regex = "^BM_OP_[SB][1-4][1-4]$", prefix = "BM_OP_"),
  BM_INV = list(regex = "^BM_INV_[SB][1-4][1-4]$", prefix = "BM_INV_"),
  OP_INV = list(regex = "^OP_INV_[SB][1-4][1-4]$", prefix = "OP_INV_"),
  OP_CHI_FIRM = list(regex = "^OP_CHI_FIRM_[SB][1-4][1-4]$", prefix = "OP_CHI_FIRM_"),
  OP_CHI_SIC = list(regex = "^OP_CHI_SIC_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC_"),
  OP_CHI_SIC2 = list(regex = "^OP_CHI_SIC2_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC2_"),
  BM_CHI_FIRM = list(regex = "^BM_CHI_FIRM_[SB][1-4][1-4]$", prefix = "BM_CHI_FIRM_"),
  BM_CHI_SIC = list(regex = "^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC_2X4X4_"),
  BM_CHI_SIC2 = list(regex = "^BM_CHI_SIC2_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC2_"),
  INV_CHI_FIRM = list(regex = "^INV_CHI_FIRM_[SB][1-4][1-4]$", prefix = "INV_CHI_FIRM_"),
  INV_CHI_SIC = list(regex = "^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC_2X4X4_"),
  INV_CHI_SIC2 = list(regex = "^INV_CHI_SIC2_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC2_")
)

pick_2x4x4_robust_models <- function(fam) {
  c(
    "FF5_SMQ_FIRM_VW",
    "FF5_SMQ_SIC_VW",
    "FF5_SMQ_SIC2_VW",
    "FF5_SMQ_FIRM_K",
    "FF5_SMQ_SIC_K",
    "FF5_SMQ_SIC2_K",
    "FF5_SMQ_FIRM_EXP_VW",
    "FF5_SMQ_SIC_EXP_VW",
    "FF5_SMQ_SIC2_EXP_VW",
    "FF5_SMQ_FIRM_EXP_K",
    "FF5_SMQ_SIC_EXP_K",
    "FF5_SMQ_SIC2_EXP_K"
  )
}

robustness_outputs_2x4x4 <- list()
robustness_alpha_regs_2x4x4 <- list()

for (fam in names(asset_blocks_2x4x4_alt)) {
  sp <- asset_blocks_2x4x4_alt[[fam]]
  out <- list()
  out_reg <- list()

  for (model_name in pick_2x4x4_robust_models(fam)) {
    reg <- run_2x4x4_regs_regex(
      factor_test_panel,
      sp$regex,
      sp$prefix,
      models_2x4x4_alt[[model_name]],
      NW_LAGS
    )

    out[[model_name]] <- list(
      S = term_to_4x4_matrix(reg, "(Intercept)", "S", digits_est = 2, digits_t = 2, scale_est = 100),
      B = term_to_4x4_matrix(reg, "(Intercept)", "B", digits_est = 2, digits_t = 2, scale_est = 100)
    )

    out_reg[[model_name]] <- extract_alpha_table(reg, portfolio_cols = c("portfolio", "size", "i", "j"))
  }

  robustness_outputs_2x4x4[[fam]] <- out
  robustness_alpha_regs_2x4x4[[fam]] <- out_reg
}

saveRDS(robustness_outputs_2x4x4, "section10_robustness_intercepts_2x4x4.rds")
saveRDS(robustness_alpha_regs_2x4x4, "section10_robustness_intercepts_2x4x4_alphas.rds")

### 10.21 Print Alternative 2x4x4 Performance #################################
display_2x4x4_tables(mean_2x4x4_alt, "SECTION 10: MONTHLY MEAN 2x4x4 RETURNS (%)")
display_2x4x4_tables(mean_excess_2x4x4_alt, "SECTION 10: MONTHLY MEAN 2x4x4 EXCESS RETURNS (%)")
display_2x4x4_tables(sd_2x4x4_alt, "SECTION 10: MONTHLY STD DEV 2x4x4 EXCESS RETURNS (%)")
display_2x4x4_tables(sr_2x4x4_alt, "SECTION 10: SHARPE RATIOS 2x4x4 (monthly)")

### 10.22 Asset Blocks Available ##############################################
asset_blocks_base <- list(
  BM_            = grep("^BM_S[1-5]V[1-5]$", names(factor_test_panel), value = TRUE),
  OP_            = grep("^OP_S[1-5]P[1-5]$", names(factor_test_panel), value = TRUE),
  INV_           = grep("^INV_S[1-5]I[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_FIRM_      = grep("^CHI_FIRM_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_FIRM_EXP_  = grep("^CHI_FIRM_EXP_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_SIC_       = grep("^CHI_SIC_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_SIC_EXP_   = grep("^CHI_SIC_EXP_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_SIC2_      = grep("^CHI_SIC2_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  CHI_SIC2_EXP_  = grep("^CHI_SIC2_EXP_S[1-5]C[1-5]$", names(factor_test_panel), value = TRUE),
  BM_OP_         = grep("^BM_OP_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  BM_INV_        = grep("^BM_INV_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  OP_INV_        = grep("^OP_INV_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  OP_CHI_FIRM_   = grep("^OP_CHI_FIRM_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  OP_CHI_SIC_    = grep("^OP_CHI_SIC_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  OP_CHI_SIC2_   = grep("^OP_CHI_SIC2_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  BM_CHI_FIRM_   = grep("^BM_CHI_FIRM_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  BM_CHI_SIC_    = grep("^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  BM_CHI_SIC2_   = grep("^BM_CHI_SIC2_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  INV_CHI_FIRM_  = grep("^INV_CHI_FIRM_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  INV_CHI_SIC_   = grep("^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE),
  INV_CHI_SIC2_  = grep("^INV_CHI_SIC2_[SB][1-4][1-4]$", names(factor_test_panel), value = TRUE)
)

### 10.23 Variant-Specific Universe Map #######################################
grs_variant_map <- list(
  SIC_VW = list(
    model = "FF5_SMQ_SIC_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_SIC_", "BM_CHI_SIC_", "INV_CHI_SIC_")
  ),
  FIRM_VW = list(
    model = "FF5_SMQ_FIRM_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_FIRM_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_FIRM_", "BM_CHI_FIRM_", "INV_CHI_FIRM_")
  ),
  SIC2_VW = list(
    model = "FF5_SMQ_SIC2_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC2_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_SIC2_", "BM_CHI_SIC2_", "INV_CHI_SIC2_")
  ),
  SIC_K = list(
    model = "FF5_SMQ_SIC_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_SIC_", "BM_CHI_SIC_", "INV_CHI_SIC_")
  ),
  FIRM_K = list(
    model = "FF5_SMQ_FIRM_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_FIRM_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_FIRM_", "BM_CHI_FIRM_", "INV_CHI_FIRM_")
  ),
  SIC2_K = list(
    model = "FF5_SMQ_SIC2_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC2_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_", "OP_CHI_SIC2_", "BM_CHI_SIC2_", "INV_CHI_SIC2_")
  ),
  SIC_EXP_VW = list(
    model = "FF5_SMQ_SIC_EXP_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  ),
  FIRM_EXP_VW = list(
    model = "FF5_SMQ_FIRM_EXP_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_FIRM_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  ),
  SIC2_EXP_VW = list(
    model = "FF5_SMQ_SIC2_EXP_VW",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC2_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  ),
  SIC_EXP_K = list(
    model = "FF5_SMQ_SIC_EXP_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  ),
  FIRM_EXP_K = list(
    model = "FF5_SMQ_FIRM_EXP_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_FIRM_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  ),
  SIC2_EXP_K = list(
    model = "FF5_SMQ_SIC2_EXP_K",
    family_25 = c("BM_", "OP_", "INV_", "CHI_SIC2_EXP_"),
    family_32 = c("BM_OP_", "BM_INV_", "OP_INV_")
  )
)

### 10.24 Optional Consistency Checks #########################################
stopifnot(length(asset_blocks_base$BM_) == 25)
stopifnot(length(asset_blocks_base$OP_) == 25)
stopifnot(length(asset_blocks_base$INV_) == 25)

stopifnot(length(asset_blocks_base$CHI_FIRM_) == 25)
stopifnot(length(asset_blocks_base$CHI_SIC_) == 25)
stopifnot(length(asset_blocks_base$CHI_SIC2_) == 25)
stopifnot(length(asset_blocks_base$CHI_FIRM_EXP_) == 25)
stopifnot(length(asset_blocks_base$CHI_SIC_EXP_) == 25)
stopifnot(length(asset_blocks_base$CHI_SIC2_EXP_) == 25)

stopifnot(length(asset_blocks_base$BM_OP_) == 32)
stopifnot(length(asset_blocks_base$BM_INV_) == 32)
stopifnot(length(asset_blocks_base$OP_INV_) == 32)

stopifnot(length(asset_blocks_base$OP_CHI_FIRM_) == 32)
stopifnot(length(asset_blocks_base$OP_CHI_SIC_) == 32)
stopifnot(length(asset_blocks_base$OP_CHI_SIC2_) == 32)
stopifnot(length(asset_blocks_base$BM_CHI_FIRM_) == 32)
stopifnot(length(asset_blocks_base$BM_CHI_SIC_) == 32)
stopifnot(length(asset_blocks_base$BM_CHI_SIC2_) == 32)
stopifnot(length(asset_blocks_base$INV_CHI_FIRM_) == 32)
stopifnot(length(asset_blocks_base$INV_CHI_SIC_) == 32)
stopifnot(length(asset_blocks_base$INV_CHI_SIC2_) == 32)

### 10.25 Expanded GRS Model Definitions ######################################
expanded_model_sets <- function(variant_name) {
  list(
    CAPM = models_general$CAPM,
    FF3 = models_general$FF3,
    FF5 = models_general$FF5,
    FF5_SMQ = variant_factor_map[[variant_name]]$ff5_smq
  )
}

### 10.26 Expanded GRS Within Matching Universes ##############################
expanded_grs_results <- list()
expanded_grs_alpha_tables <- list()

for (variant_name in names(grs_variant_map)) {
  spec <- grs_variant_map[[variant_name]]
  model_map <- expanded_model_sets(variant_name)

  for (model_label in names(model_map)) {
    factor_cols <- model_map[[model_label]]

    for (block_name in spec$family_25) {
      asset_cols <- asset_blocks_base[[block_name]]
      if (length(asset_cols) == 0) next

      block_res <- run_block_model_with_alphas(
        dt = factor_test_panel,
        asset_cols = asset_cols,
        factor_cols = factor_cols,
        lags = NW_LAGS
      )

      if (!is.null(block_res$summary)) {
        res <- copy(block_res$summary)
        res[, variant := variant_name]
        res[, asset_group := "family_25"]
        res[, asset_block := block_name]
        res[, model := model_label]
        expanded_grs_results[[paste(variant_name, model_label, "family25", block_name, sep = "_")]] <- res
      }

      alpha_dt <- copy(block_res$alphas)
      alpha_dt[, variant := variant_name]
      alpha_dt[, asset_group := "family_25"]
      alpha_dt[, asset_block := block_name]
      alpha_dt[, model := model_label]
      expanded_grs_alpha_tables[[paste(variant_name, model_label, "family25", block_name, sep = "_")]] <- alpha_dt
    }

    for (block_name in spec$family_32) {
      asset_cols <- asset_blocks_base[[block_name]]
      if (length(asset_cols) == 0) next

      block_res <- run_block_model_with_alphas(
        dt = factor_test_panel,
        asset_cols = asset_cols,
        factor_cols = factor_cols,
        lags = NW_LAGS
      )

      if (!is.null(block_res$summary)) {
        res <- copy(block_res$summary)
        res[, variant := variant_name]
        res[, asset_group := "family_32"]
        res[, asset_block := block_name]
        res[, model := model_label]
        expanded_grs_results[[paste(variant_name, model_label, "family32", block_name, sep = "_")]] <- res
      }

      alpha_dt <- copy(block_res$alphas)
      alpha_dt[, variant := variant_name]
      alpha_dt[, asset_group := "family_32"]
      alpha_dt[, asset_block := block_name]
      alpha_dt[, model := model_label]
      expanded_grs_alpha_tables[[paste(variant_name, model_label, "family32", block_name, sep = "_")]] <- alpha_dt
    }

    all_25_cols <- unique(unlist(asset_blocks_base[spec$family_25]))
    if (length(all_25_cols) > 0) {
      block_res <- run_block_model_with_alphas(
        dt = factor_test_panel,
        asset_cols = all_25_cols,
        factor_cols = factor_cols,
        lags = NW_LAGS
      )

      if (!is.null(block_res$summary)) {
        res <- copy(block_res$summary)
        res[, variant := variant_name]
        res[, asset_group := "all_25"]
        res[, asset_block := "ALL_25"]
        res[, model := model_label]
        expanded_grs_results[[paste(variant_name, model_label, "ALL_25", sep = "_")]] <- res
      }

      alpha_dt <- copy(block_res$alphas)
      alpha_dt[, variant := variant_name]
      alpha_dt[, asset_group := "all_25"]
      alpha_dt[, asset_block := "ALL_25"]
      alpha_dt[, model := model_label]
      expanded_grs_alpha_tables[[paste(variant_name, model_label, "ALL_25", sep = "_")]] <- alpha_dt
    }

    all_32_cols <- unique(unlist(asset_blocks_base[spec$family_32]))
    if (length(all_32_cols) > 0) {
      block_res <- run_block_model_with_alphas(
        dt = factor_test_panel,
        asset_cols = all_32_cols,
        factor_cols = factor_cols,
        lags = NW_LAGS
      )

      if (!is.null(block_res$summary)) {
        res <- copy(block_res$summary)
        res[, variant := variant_name]
        res[, asset_group := "all_32"]
        res[, asset_block := "ALL_32"]
        res[, model := model_label]
        expanded_grs_results[[paste(variant_name, model_label, "ALL_32", sep = "_")]] <- res
      }

      alpha_dt <- copy(block_res$alphas)
      alpha_dt[, variant := variant_name]
      alpha_dt[, asset_group := "all_32"]
      alpha_dt[, asset_block := "ALL_32"]
      alpha_dt[, model := model_label]
      expanded_grs_alpha_tables[[paste(variant_name, model_label, "ALL_32", sep = "_")]] <- alpha_dt
    }

    all_variant_cols <- unique(c(all_25_cols, all_32_cols))
    if (length(all_variant_cols) > 0) {
      block_res <- run_block_model_with_alphas(
        dt = factor_test_panel,
        asset_cols = all_variant_cols,
        factor_cols = factor_cols,
        lags = NW_LAGS
      )

      if (!is.null(block_res$summary)) {
        res <- copy(block_res$summary)
        res[, variant := variant_name]
        res[, asset_group := "all_assets"]
        res[, asset_block := "ALL_ASSETS"]
        res[, model := model_label]
        expanded_grs_results[[paste(variant_name, model_label, "ALL_ASSETS", sep = "_")]] <- res
      }

      alpha_dt <- copy(block_res$alphas)
      alpha_dt[, variant := variant_name]
      alpha_dt[, asset_group := "all_assets"]
      alpha_dt[, asset_block := "ALL_ASSETS"]
      alpha_dt[, model := model_label]
      expanded_grs_alpha_tables[[paste(variant_name, model_label, "ALL_ASSETS", sep = "_")]] <- alpha_dt
    }
  }
}

expanded_grs_summary_table <- rbindlist(expanded_grs_results, fill = TRUE)
setcolorder(
  expanded_grs_summary_table,
  c(
    "variant", "asset_group", "asset_block", "model",
    "T", "N", "K", "GRS", "GRS_pval",
    "mean_abs_alpha", "rms_alpha", "avg_adj_r2",
    "A_abs_a_over_A_abs_rbar", "A_a2_over_A_rbar2"
  )
)

expanded_grs_alpha_table <- rbindlist(expanded_grs_alpha_tables, fill = TRUE)
saveRDS(expanded_grs_summary_table, "section10_expanded_grs_summary_table.rds")
saveRDS(expanded_grs_alpha_table, "section10_expanded_grs_alpha_table.rds")

### 10.27 Average Adjusted R2 Summary #########################################
avg_adj_r2_summary <- expanded_grs_summary_table[
  ,
  .(
    avg_adj_r2 = mean(avg_adj_r2, na.rm = TRUE),
    avg_mean_abs_alpha = mean(mean_abs_alpha, na.rm = TRUE),
    avg_rms_alpha = mean(rms_alpha, na.rm = TRUE),
    avg_GRS = mean(GRS, na.rm = TRUE),
    avg_GRS_pval = mean(GRS_pval, na.rm = TRUE)
  ),
  by = .(variant, model, asset_group)
][order(variant, model, asset_group)]

saveRDS(avg_adj_r2_summary, "section10_avg_adj_r2_summary_table.rds")

### 10.28 CAPM-Only Chi Portfolio Alpha Outputs ###############################
capm_only_chi_alpha_tables <- list()
capm_only_chi_significance <- list()
capm_vs_other_alpha_count_comparison <- list()

run_capm_alpha_block_5x5 <- function(regex_pattern, strip_prefix) {
  reg <- run_asset_block_regs_regex(
    panel = factor_test_panel,
    regex_pattern = regex_pattern,
    strip_prefix = strip_prefix,
    rhs = models_general$CAPM,
    lags = NW_LAGS
  )
  extract_alpha_table(reg)
}

run_model_alpha_block_5x5 <- function(regex_pattern, strip_prefix, rhs) {
  reg <- run_asset_block_regs_regex(
    panel = factor_test_panel,
    regex_pattern = regex_pattern,
    strip_prefix = strip_prefix,
    rhs = rhs,
    lags = NW_LAGS
  )
  extract_alpha_table(reg)
}

run_capm_alpha_block_2x4x4 <- function(regex_pattern, strip_prefix) {
  reg <- run_2x4x4_regs_regex(
    panel = factor_test_panel,
    regex_pattern = regex_pattern,
    strip_prefix = strip_prefix,
    rhs = models_general$CAPM,
    lags = NW_LAGS
  )
  extract_alpha_table(reg, portfolio_cols = c("portfolio", "size", "i", "j"))
}

run_model_alpha_block_2x4x4 <- function(regex_pattern, strip_prefix, rhs) {
  reg <- run_2x4x4_regs_regex(
    panel = factor_test_panel,
    regex_pattern = regex_pattern,
    strip_prefix = strip_prefix,
    rhs = rhs,
    lags = NW_LAGS
  )
  extract_alpha_table(reg, portfolio_cols = c("portfolio", "size", "i", "j"))
}

chi_25_specs <- list(
  CHI_FIRM_25 = list(regex = "^CHI_FIRM_S[1-5]C[1-5]$", prefix = "CHI_FIRM_"),
  CHI_SIC_25 = list(regex = "^CHI_SIC_S[1-5]C[1-5]$", prefix = "CHI_SIC_"),
  CHI_SIC2_25 = list(regex = "^CHI_SIC2_S[1-5]C[1-5]$", prefix = "CHI_SIC2_"),
  CHI_FIRM_EXP_25 = list(regex = "^CHI_FIRM_EXP_S[1-5]C[1-5]$", prefix = "CHI_FIRM_EXP_"),
  CHI_SIC_EXP_25 = list(regex = "^CHI_SIC_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC_EXP_"),
  CHI_SIC2_EXP_25 = list(regex = "^CHI_SIC2_EXP_S[1-5]C[1-5]$", prefix = "CHI_SIC2_EXP_")
)

for (fam in names(chi_25_specs)) {
  sp <- chi_25_specs[[fam]]

  capm_alpha <- run_capm_alpha_block_5x5(sp$regex, sp$prefix)
  capm_alpha[, family := fam]
  capm_alpha[, chi_bucket := char]
  capm_alpha[, size_bucket := size]

  capm_only_chi_alpha_tables[[fam]] <- capm_alpha

  capm_only_chi_significance[[fam]] <- summarize_significant_alphas(
    alpha_dt = capm_alpha,
    family_name = fam,
    chi_pos = "chi_bucket",
    size_pos = "size_bucket"
  )
}

chi_32_specs <- list(
  OP_CHI_FIRM  = list(regex = "^OP_CHI_FIRM_[SB][1-4][1-4]$", prefix = "OP_CHI_FIRM_"),
  OP_CHI_SIC   = list(regex = "^OP_CHI_SIC_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC_"),
  OP_CHI_SIC2  = list(regex = "^OP_CHI_SIC2_[SB][1-4][1-4]$", prefix = "OP_CHI_SIC2_"),
  BM_CHI_FIRM  = list(regex = "^BM_CHI_FIRM_[SB][1-4][1-4]$", prefix = "BM_CHI_FIRM_"),
  BM_CHI_SIC   = list(regex = "^BM_CHI_SIC_2X4X4_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC_2X4X4_"),
  BM_CHI_SIC2  = list(regex = "^BM_CHI_SIC2_[SB][1-4][1-4]$", prefix = "BM_CHI_SIC2_"),
  INV_CHI_FIRM = list(regex = "^INV_CHI_FIRM_[SB][1-4][1-4]$", prefix = "INV_CHI_FIRM_"),
  INV_CHI_SIC  = list(regex = "^INV_CHI_SIC_2X4X4_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC_2X4X4_"),
  INV_CHI_SIC2 = list(regex = "^INV_CHI_SIC2_[SB][1-4][1-4]$", prefix = "INV_CHI_SIC2_")
)

for (fam in names(chi_32_specs)) {
  sp <- chi_32_specs[[fam]]

  capm_alpha <- run_capm_alpha_block_2x4x4(sp$regex, sp$prefix)
  capm_alpha[, family := fam]
  capm_alpha[, chi_bucket := j]
  capm_alpha[, size_bucket := size]

  capm_only_chi_alpha_tables[[fam]] <- capm_alpha

  capm_only_chi_significance[[fam]] <- summarize_significant_alphas(
    alpha_dt = capm_alpha,
    family_name = fam,
    chi_pos = "chi_bucket",
    size_pos = "size_bucket"
  )
}

capm_only_chi_alpha_table <- rbindlist(capm_only_chi_alpha_tables, fill = TRUE)
capm_only_chi_significance_summary <- rbindlist(capm_only_chi_significance, fill = TRUE)

saveRDS(capm_only_chi_alpha_tables, "section10_capm_chi_alpha_tables_list.rds")
saveRDS(capm_only_chi_alpha_table, "section10_capm_chi_alpha_table.rds")
saveRDS(capm_only_chi_significance_summary, "section10_capm_chi_significance_summary.rds")

### 10.29 CAPM vs FF3 vs FF5 vs FF5+SMQ Alpha Count Metrics ###################
chi_family_model_counts <- list()

chi_family_variant_for_smq <- function(fam) {
  if (grepl("FIRM", fam)) {
    return("FIRM_VW")
  }
  if (grepl("SIC2", fam)) {
    return("SIC2_VW")
  }
  if (grepl("SIC", fam)) {
    return("SIC_VW")
  }
  NA_character_
}

for (fam in names(chi_32_specs)) {
  sp <- chi_32_specs[[fam]]
  variant_name <- chi_family_variant_for_smq(fam)
  ff5_smq_rhs <- variant_factor_map[[variant_name]]$ff5_smq

  alpha_list <- list(
    CAPM = {
      tmp <- run_capm_alpha_block_2x4x4(sp$regex, sp$prefix)
      tmp[, chi_bucket := j]
      tmp[, size_bucket := size]
      tmp
    },
    FF3 = {
      tmp <- run_model_alpha_block_2x4x4(sp$regex, sp$prefix, models_general$FF3)
      tmp[, chi_bucket := j]
      tmp[, size_bucket := size]
      tmp
    },
    FF5 = {
      tmp <- run_model_alpha_block_2x4x4(sp$regex, sp$prefix, models_general$FF5)
      tmp[, chi_bucket := j]
      tmp[, size_bucket := size]
      tmp
    },
    FF5_SMQ = {
      tmp <- run_model_alpha_block_2x4x4(sp$regex, sp$prefix, ff5_smq_rhs)
      tmp[, chi_bucket := j]
      tmp[, size_bucket := size]
      tmp
    }
  )

  chi_family_model_counts[[fam]] <- make_alpha_count_comparison(
    alpha_list_named = alpha_list,
    family_name = fam,
    size_col = "size_bucket",
    chi_col = "chi_bucket"
  )
}

chi_25_model_counts <- list()

for (fam in names(chi_25_specs)) {
  sp <- chi_25_specs[[fam]]
  variant_name <- chi_family_variant_for_smq(fam)
  ff5_smq_rhs <- variant_factor_map[[variant_name]]$ff5_smq

  alpha_list <- list(
    CAPM = {
      tmp <- run_capm_alpha_block_5x5(sp$regex, sp$prefix)
      tmp[, chi_bucket := char]
      tmp[, size_bucket := size]
      tmp
    },
    FF3 = {
      tmp <- run_model_alpha_block_5x5(sp$regex, sp$prefix, models_general$FF3)
      tmp[, chi_bucket := char]
      tmp[, size_bucket := size]
      tmp
    },
    FF5 = {
      tmp <- run_model_alpha_block_5x5(sp$regex, sp$prefix, models_general$FF5)
      tmp[, chi_bucket := char]
      tmp[, size_bucket := size]
      tmp
    },
    FF5_SMQ = {
      tmp <- run_model_alpha_block_5x5(sp$regex, sp$prefix, ff5_smq_rhs)
      tmp[, chi_bucket := char]
      tmp[, size_bucket := size]
      tmp
    }
  )

  chi_25_model_counts[[fam]] <- make_alpha_count_comparison(
    alpha_list_named = alpha_list,
    family_name = fam,
    size_col = "size_bucket",
    chi_col = "chi_bucket"
  )
}

chi_family_model_counts_table <- rbindlist(c(chi_family_model_counts, chi_25_model_counts), fill = TRUE)
saveRDS(chi_family_model_counts_table, "section10_chi_family_model_count_comparison.rds")

### 10.30 Alternative Subperiod Evidence ######################################
subperiods <- list(
  Early = list(start = as.Date("1982-07-01"), end = as.Date("2002-12-01")),
  Late  = list(start = as.Date("2003-01-01"), end = as.Date("2024-12-01"))
)

alt_subperiod_factors <- smq_variants

return_panel_alt <- merge(
  ff5[, .(mdate, HML, CMA)],
  SMQ[, c("mdate", alt_subperiod_factors), with = FALSE],
  by = "mdate",
  all = FALSE
)

return_panel_alt[, period := fifelse(
  mdate <= as.Date("2002-12-01"),
  "Early",
  "Late"
)]

alt_subperiod_stats <- list()
alt_subperiod_corrs <- list()
alt_subperiod_factor_regs <- list()
alt_rolling_panels <- list()

for (fac in alt_subperiod_factors) {
  alt_subperiod_stats[[fac]] <- list()
  alt_subperiod_corrs[[fac]] <- list()
  alt_subperiod_factor_regs[[fac]] <- list()

  for (nm in names(subperiods)) {
    sp <- subperiods[[nm]]

    dt_p <- return_panel_alt[
      mdate >= sp$start &
        mdate <= sp$end,
      .(mdate, HML, CMA, SMQ = get(fac))
    ]

    alt_subperiod_stats[[fac]][[nm]] <- make_subperiod_summary_table(
      dt_p,
      c("HML", "CMA", "SMQ"),
      NW_LAGS
    )

    alt_subperiod_corrs[[fac]][[nm]] <- cor(
      dt_p[, .(HML, CMA, SMQ)],
      use = "pairwise.complete.obs"
    )

    alt_subperiod_factor_regs[[fac]][[nm]] <- make_subperiod_factor_regression_table(
      dt_p,
      smq_col = "SMQ",
      lags = NW_LAGS
    )
  }

  plot_panel <- return_panel_alt[, .(mdate, HML, CMA, fac_value = get(fac))]
  plot_panel[, HML_roll := roll_mean_safe(HML, ROLL_WIN)]
  plot_panel[, CMA_roll := roll_mean_safe(CMA, ROLL_WIN)]
  plot_panel[, FAC_roll := roll_mean_safe(fac_value, ROLL_WIN)]
  plot_panel[, HML_sr_roll := roll_sharpe_safe(HML, ROLL_WIN)]
  plot_panel[, CMA_sr_roll := roll_sharpe_safe(CMA, ROLL_WIN)]
  plot_panel[, FAC_sr_roll := roll_sharpe_safe(fac_value, ROLL_WIN)]
  alt_rolling_panels[[fac]] <- plot_panel
}

alt_subperiod_stats_table <- rbindlist(
  lapply(names(alt_subperiod_stats), function(fac) {
    rbindlist(
      lapply(names(alt_subperiod_stats[[fac]]), function(prd) {
        tmp <- copy(alt_subperiod_stats[[fac]][[prd]])
        tmp[, factor := fac]
        tmp[, period := prd]
        tmp
      }),
      fill = TRUE
    )
  }),
  fill = TRUE
)

saveRDS(alt_subperiod_stats, "section10_alt_subperiod_stats_list.rds")
saveRDS(alt_subperiod_stats_table, "section10_alt_subperiod_stats.rds")
saveRDS(alt_subperiod_corrs, "section10_alt_subperiod_correlations.rds")
saveRDS(alt_subperiod_factor_regs, "section10_alt_subperiod_factor_regs.rds")
saveRDS(alt_rolling_panels, "section10_alt_rolling_panels.rds")

### 10.31 Alternative Rolling Plots ###########################################
for (fac in alt_subperiod_factors) {
  plot_panel <- alt_rolling_panels[[fac]]

  par(mfrow = c(1, 1))
  plot(
    plot_panel$mdate,
    plot_panel$HML_roll,
    type = "l",
    lwd = 2,
    main = paste("60-Month Rolling Means:", fac),
    xlab = "",
    ylab = "Rolling Mean Return"
  )
  lines(plot_panel$mdate, plot_panel$CMA_roll, lwd = 2, lty = 2)
  lines(plot_panel$mdate, plot_panel$FAC_roll, lwd = 2, lty = 3)
  abline(v = as.Date("2003-01-01"), lty = 4)
  legend(
    "topright",
    legend = c("HML", "CMA", fac, "Early/Late split"),
    lty = c(1, 2, 3, 4),
    lwd = c(2, 2, 2, 1),
    bty = "n"
  )

  par(mfrow = c(1, 1))
  plot(
    plot_panel$mdate,
    plot_panel$HML_sr_roll,
    type = "l",
    lwd = 2,
    main = paste("60-Month Rolling Sharpe Ratios:", fac),
    xlab = "",
    ylab = "Rolling Sharpe Ratio"
  )
  lines(plot_panel$mdate, plot_panel$CMA_sr_roll, lwd = 2, lty = 2)
  lines(plot_panel$mdate, plot_panel$FAC_sr_roll, lwd = 2, lty = 3)
  abline(v = as.Date("2003-01-01"), lty = 4)
  legend(
    "topright",
    legend = c("HML", "CMA", fac, "Early/Late split"),
    lty = c(1, 2, 3, 4),
    lwd = c(2, 2, 2, 1),
    bty = "n"
  )
}

### 10.32 Additional Explicit CAPM Display Objects ############################
make_portfolio_alpha_display <- function(alpha_dt, scale = 100) {
  out <- copy(alpha_dt)
  out[, alpha_display := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    scale * alpha, alpha_stars, alpha_tstat, alpha_pval
  )]
  out[]
}

capm_chi_alpha_display_list <- lapply(capm_only_chi_alpha_tables, make_portfolio_alpha_display)
saveRDS(capm_chi_alpha_display_list, "section10_capm_chi_alpha_display_list.rds")

### 10.33 Expanded GRS Pivot-Friendly Summary #################################
expanded_grs_pivot_summary <- expanded_grs_summary_table[
  ,
  .(
    variant,
    asset_group,
    asset_block,
    model,
    GRS,
    GRS_pval,
    mean_abs_alpha,
    rms_alpha,
    avg_adj_r2,
    A_abs_a_over_A_abs_rbar,
    A_a2_over_A_rbar2
  )
][order(variant, asset_group, asset_block, model)]

saveRDS(expanded_grs_pivot_summary, "section10_expanded_grs_pivot_summary.rds")

### 10.34 Model-by-Variant Comparison Tables ##################################
model_variant_comparison_table <- expanded_grs_summary_table[
  ,
  .(
    avg_GRS = mean(GRS, na.rm = TRUE),
    med_GRS = median(GRS, na.rm = TRUE),
    avg_GRS_pval = mean(GRS_pval, na.rm = TRUE),
    avg_mean_abs_alpha = mean(mean_abs_alpha, na.rm = TRUE),
    avg_rms_alpha = mean(rms_alpha, na.rm = TRUE),
    avg_adj_r2 = mean(avg_adj_r2, na.rm = TRUE)
  ),
  by = .(variant, model)
][order(variant, model)]

saveRDS(model_variant_comparison_table, "section10_model_variant_comparison_table.rds")

### 10.35 CAPM Chi Family Summary Tables ######################################
capm_chi_family_summary <- capm_only_chi_alpha_table[
  ,
  .(
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_abs_alpha = mean(abs(alpha), na.rm = TRUE),
    median_abs_alpha = median(abs(alpha), na.rm = TRUE),
    share_sig_10 = mean(alpha_pval < 0.10, na.rm = TRUE),
    share_sig_05 = mean(alpha_pval < 0.05, na.rm = TRUE),
    share_sig_01 = mean(alpha_pval < 0.01, na.rm = TRUE),
    n_pos_sig_05 = sum(alpha_pval < 0.05 & alpha > 0, na.rm = TRUE),
    n_neg_sig_05 = sum(alpha_pval < 0.05 & alpha < 0, na.rm = TRUE)
  ),
  by = .(family)
][order(family)]

saveRDS(capm_chi_family_summary, "section10_capm_chi_family_summary.rds")

### 10.36 CAPM Chi Bucket Tables ##############################################
capm_chi_bucket_summary <- capm_only_chi_alpha_table[
  ,
  .(
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_abs_alpha = mean(abs(alpha), na.rm = TRUE),
    share_sig_10 = mean(alpha_pval < 0.10, na.rm = TRUE),
    share_sig_05 = mean(alpha_pval < 0.05, na.rm = TRUE),
    share_sig_01 = mean(alpha_pval < 0.01, na.rm = TRUE),
    n_pos_sig_05 = sum(alpha_pval < 0.05 & alpha > 0, na.rm = TRUE),
    n_neg_sig_05 = sum(alpha_pval < 0.05 & alpha < 0, na.rm = TRUE)
  ),
  by = .(family, size_bucket, chi_bucket)
][order(family, size_bucket, chi_bucket)]

saveRDS(capm_chi_bucket_summary, "section10_capm_chi_bucket_summary.rds")

### 10.37 Final Saved Master Objects ##########################################
section10_master_outputs <- list(
  smq_summary_table = smq_summary_table,
  smq_summary_display = smq_summary_display,
  expanded_grs_summary_table = expanded_grs_summary_table,
  expanded_grs_alpha_table = expanded_grs_alpha_table,
  avg_adj_r2_summary = avg_adj_r2_summary,
  capm_only_chi_alpha_table = capm_only_chi_alpha_table,
  capm_only_chi_significance_summary = capm_only_chi_significance_summary,
  chi_family_model_counts_table = chi_family_model_counts_table,
  model_variant_comparison_table = model_variant_comparison_table,
  capm_chi_family_summary = capm_chi_family_summary,
  capm_chi_bucket_summary = capm_chi_bucket_summary
)

saveRDS(section10_master_outputs, "section10_master_outputs.rds")

### 10.38 Console Prints ######################################################
cat("\n============================================================\n")
cat("SECTION 10 ROBUSTNESS COMPLETE\n")
cat("============================================================\n\n")

cat("Saved key outputs:\n")
cat("- section10_smq_summary_table.rds\n")
cat("- section10_expanded_grs_summary_table.rds\n")
cat("- section10_expanded_grs_alpha_table.rds\n")
cat("- section10_avg_adj_r2_summary_table.rds\n")
cat("- section10_capm_chi_alpha_table.rds\n")
cat("- section10_capm_chi_significance_summary.rds\n")
cat("- section10_chi_family_model_count_comparison.rds\n")
cat("- section10_master_outputs.rds\n\n")

print(round_dt(smq_summary_table, 4))
print(round_dt(expanded_grs_summary_table, 4))
print(round_dt(avg_adj_r2_summary, 4))
print(round_dt(capm_chi_family_summary, 4))






### 10.39 Combined Requested Output and Multi-Model Display Section ###########

### 10.39A Setup ##############################################################
capm_rhs <- c("MKT")
ff3_rhs <- c("MKT", "SMB", "HML")
ff5_rhs <- c("MKT", "SMB", "HML", "RMW", "CMA")

requested_models <- c("CAPM", "FF3", "FF5", "FF5_SMQ")
model_order <- requested_models
variant_order <- c("FIRM_VW", "SIC_VW", "SIC2_VW")

chi_25_order <- c("CHI_FIRM_25", "CHI_SIC_25", "CHI_SIC2_25")

chi_32_order <- c(
  "OP_CHI_FIRM_32", "BM_CHI_FIRM_32", "INV_CHI_FIRM_32",
  "OP_CHI_SIC_32",  "BM_CHI_SIC_32",  "INV_CHI_SIC_32",
  "OP_CHI_SIC2_32", "BM_CHI_SIC2_32", "INV_CHI_SIC2_32"
)

chi_25_block_order <- c("CHI_FIRM_", "CHI_SIC_", "CHI_SIC2_")

chi_32_block_order <- c(
  "OP_CHI_FIRM_", "BM_CHI_FIRM_", "INV_CHI_FIRM_",
  "OP_CHI_SIC_",  "BM_CHI_SIC_",  "INV_CHI_SIC_",
  "OP_CHI_SIC2_", "BM_CHI_SIC2_", "INV_CHI_SIC2_"
)

chi_variant_lookup_25 <- data.table(
  asset_block = c("CHI_FIRM_", "CHI_SIC_", "CHI_SIC2_"),
  variant = c("FIRM_VW", "SIC_VW", "SIC2_VW"),
  portfolio_set = c("CHI_FIRM_25", "CHI_SIC_25", "CHI_SIC2_25")
)

chi_variant_lookup_32 <- data.table(
  asset_block = c(
    "OP_CHI_FIRM_", "BM_CHI_FIRM_", "INV_CHI_FIRM_",
    "OP_CHI_SIC_",  "BM_CHI_SIC_",  "INV_CHI_SIC_",
    "OP_CHI_SIC2_", "BM_CHI_SIC2_", "INV_CHI_SIC2_"
  ),
  variant = c(
    "FIRM_VW", "FIRM_VW", "FIRM_VW",
    "SIC_VW",  "SIC_VW",  "SIC_VW",
    "SIC2_VW", "SIC2_VW", "SIC2_VW"
  ),
  portfolio_set = c(
    "OP_CHI_FIRM_32", "BM_CHI_FIRM_32", "INV_CHI_FIRM_32",
    "OP_CHI_SIC_32",  "BM_CHI_SIC_32",  "INV_CHI_SIC_32",
    "OP_CHI_SIC2_32", "BM_CHI_SIC2_32", "INV_CHI_SIC2_32"
  )
)

### 10.39A1 Robust portfolio parsers ##########################################
parse_chi_25_portfolio <- function(portfolio) {
  out <- data.table(
    portfolio = portfolio,
    size_bucket = NA_character_,
    chi_bucket = NA_character_
  )

  x <- sub("^CHI_(FIRM|SIC|SIC2)(_EXP)?_", "", portfolio)
  m <- regexec("^(S[1-5])C([1-5])$", x)
  g <- regmatches(x, m)[[1]]

  if (length(g) == 3) {
    out[, size_bucket := g[2]]
    out[, chi_bucket := paste0("C", g[3])]
  }

  out
}

parse_chi_32_portfolio <- function(portfolio) {
  out <- data.table(
    portfolio = portfolio,
    size_bucket = NA_character_,
    char_bucket = NA_character_,
    chi_bucket = NA_character_
  )

  x <- sub("^(OP|BM|INV)_CHI_(FIRM|SIC|SIC2)(_2X4X4)?_", "", portfolio)
  m <- regexec("^([SB])(\\d)(\\d)$", x)
  g <- regmatches(x, m)[[1]]

  if (length(g) == 4) {
    out[, size_bucket := g[2]]
    out[, char_bucket := paste0("Q", g[3])]
    out[, chi_bucket := paste0("Q", g[4])]
  }

  out
}

attach_chi_25_labels <- function(dt) {
  parsed <- rbindlist(
    lapply(unique(dt$portfolio), parse_chi_25_portfolio),
    fill = TRUE
  )
  parsed <- parsed[!duplicated(portfolio)]
  merge(dt, parsed, by = "portfolio", all.x = TRUE)
}

attach_chi_32_labels <- function(dt) {
  parsed <- rbindlist(
    lapply(unique(dt$portfolio), parse_chi_32_portfolio),
    fill = TRUE
  )
  parsed <- parsed[!duplicated(portfolio)]
  merge(dt, parsed, by = "portfolio", all.x = TRUE)
}

### 10.39B Multi-model GRS displays ###########################################

## 25-portfolio chi GRS blocks
model_grs_chi_25_vw <- merge(
  expanded_grs_summary_table[
    model %in% requested_models &
      asset_group == "family_25" &
      asset_block %in% chi_25_block_order &
      variant %in% variant_order
  ],
  chi_variant_lookup_25,
  by = c("asset_block", "variant"),
  all.x = TRUE
)

model_grs_chi_25_vw <- model_grs_chi_25_vw[
  order(
    match(portfolio_set, chi_25_order),
    match(model, model_order)
  )
]

model_grs_chi_25_vw_display <- copy(model_grs_chi_25_vw)

setcolorder(
  model_grs_chi_25_vw_display,
  c(
    "variant", "model", "asset_group", "asset_block", "portfolio_set",
    "T", "N", "K", "GRS", "GRS_pval",
    "mean_abs_alpha", "rms_alpha", "avg_adj_r2",
    "A_abs_a_over_A_abs_rbar", "A_a2_over_A_rbar2"
  )
)

saveRDS(model_grs_chi_25_vw, "section10_model_grs_chi_25_vw.rds")
saveRDS(model_grs_chi_25_vw_display, "section10_model_grs_chi_25_vw_display.rds")

## 32-portfolio chi GRS blocks
model_grs_chi_32_vw <- merge(
  expanded_grs_summary_table[
    model %in% requested_models &
      asset_group == "family_32" &
      asset_block %in% chi_32_block_order &
      variant %in% variant_order
  ],
  chi_variant_lookup_32,
  by = c("asset_block", "variant"),
  all.x = TRUE
)

model_grs_chi_32_vw <- model_grs_chi_32_vw[
  order(
    match(portfolio_set, chi_32_order),
    match(model, model_order)
  )
]

model_grs_chi_32_vw_display <- copy(model_grs_chi_32_vw)

setcolorder(
  model_grs_chi_32_vw_display,
  c(
    "variant", "model", "asset_group", "asset_block", "portfolio_set",
    "T", "N", "K", "GRS", "GRS_pval",
    "mean_abs_alpha", "rms_alpha", "avg_adj_r2",
    "A_abs_a_over_A_abs_rbar", "A_a2_over_A_rbar2"
  )
)

saveRDS(model_grs_chi_32_vw, "section10_model_grs_chi_32_vw.rds")
saveRDS(model_grs_chi_32_vw_display, "section10_model_grs_chi_32_vw_display.rds")

### 10.39C Multi-model individual alpha displays ##############################

## 25-portfolio chi individual alpha table
model_alpha_chi_25_vw <- merge(
  expanded_grs_alpha_table[
    model %in% requested_models &
      asset_group == "family_25" &
      asset_block %in% chi_25_block_order &
      variant %in% variant_order
  ],
  chi_variant_lookup_25,
  by = c("asset_block", "variant"),
  all.x = TRUE
)

model_alpha_chi_25_vw <- attach_chi_25_labels(model_alpha_chi_25_vw)
model_alpha_chi_25_vw[, alpha_pct := 100 * alpha]
model_alpha_chi_25_vw[, alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  alpha_pct, alpha_stars, alpha_tstat, alpha_pval
)]

model_alpha_chi_25_vw <- model_alpha_chi_25_vw[
  order(
    match(portfolio_set, chi_25_order),
    match(model, model_order),
    size_bucket,
    chi_bucket
  )
]

model_alpha_chi_25_vw_display <- copy(model_alpha_chi_25_vw)

setcolorder(
  model_alpha_chi_25_vw_display,
  c(
    "variant", "model", "portfolio_set", "portfolio",
    "size_bucket", "chi_bucket",
    "alpha", "alpha_pct", "alpha_tstat", "alpha_pval", "alpha_stars", "adj_r2",
    "alpha_display"
  )
)

saveRDS(model_alpha_chi_25_vw, "section10_model_alpha_chi_25_vw.rds")
saveRDS(model_alpha_chi_25_vw_display, "section10_model_alpha_chi_25_vw_display.rds")

## 32-portfolio chi individual alpha table
model_alpha_chi_32_vw <- merge(
  expanded_grs_alpha_table[
    model %in% requested_models &
      asset_group == "family_32" &
      asset_block %in% chi_32_block_order &
      variant %in% variant_order
  ],
  chi_variant_lookup_32,
  by = c("asset_block", "variant"),
  all.x = TRUE
)

model_alpha_chi_32_vw <- attach_chi_32_labels(model_alpha_chi_32_vw)
model_alpha_chi_32_vw[, alpha_pct := 100 * alpha]
model_alpha_chi_32_vw[, alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  alpha_pct, alpha_stars, alpha_tstat, alpha_pval
)]

model_alpha_chi_32_vw <- model_alpha_chi_32_vw[
  order(
    match(portfolio_set, chi_32_order),
    match(model, model_order),
    size_bucket,
    char_bucket,
    chi_bucket
  )
]

model_alpha_chi_32_vw_display <- copy(model_alpha_chi_32_vw)

setcolorder(
  model_alpha_chi_32_vw_display,
  c(
    "variant", "model", "portfolio_set", "portfolio",
    "size_bucket", "char_bucket", "chi_bucket",
    "alpha", "alpha_pct", "alpha_tstat", "alpha_pval", "alpha_stars", "adj_r2",
    "alpha_display"
  )
)

saveRDS(model_alpha_chi_32_vw, "section10_model_alpha_chi_32_vw.rds")
saveRDS(model_alpha_chi_32_vw_display, "section10_model_alpha_chi_32_vw_display.rds")

### 10.39D Multi-model alpha summaries ########################################
## 25-portfolio chi alpha summary
model_alpha_chi_25_vw_summary <- model_alpha_chi_25_vw[
  ,
  .(
    n_portfolios = .N,
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_abs_alpha = mean(abs(alpha), na.rm = TRUE),
    median_abs_alpha = median(abs(alpha), na.rm = TRUE),
    share_sig_10 = mean(alpha_pval < 0.10, na.rm = TRUE),
    share_sig_05 = mean(alpha_pval < 0.05, na.rm = TRUE),
    share_sig_01 = mean(alpha_pval < 0.01, na.rm = TRUE),
    n_sig_05 = sum(alpha_pval < 0.05, na.rm = TRUE),
    n_pos_sig_05 = sum(alpha_pval < 0.05 & alpha > 0, na.rm = TRUE),
    n_neg_sig_05 = sum(alpha_pval < 0.05 & alpha < 0, na.rm = TRUE),
    avg_adj_r2 = mean(adj_r2, na.rm = TRUE)
  ),
  by = .(portfolio_set, variant, model)
][order(
  match(portfolio_set, chi_25_order),
  match(model, model_order)
)]

saveRDS(
  model_alpha_chi_25_vw_summary,
  "section10_model_alpha_chi_25_vw_summary.rds"
)

## 32-portfolio chi alpha summary
model_alpha_chi_32_vw_summary <- model_alpha_chi_32_vw[
  ,
  .(
    n_portfolios = .N,
    mean_alpha = mean(alpha, na.rm = TRUE),
    mean_abs_alpha = mean(abs(alpha), na.rm = TRUE),
    median_abs_alpha = median(abs(alpha), na.rm = TRUE),
    share_sig_10 = mean(alpha_pval < 0.10, na.rm = TRUE),
    share_sig_05 = mean(alpha_pval < 0.05, na.rm = TRUE),
    share_sig_01 = mean(alpha_pval < 0.01, na.rm = TRUE),
    n_sig_05 = sum(alpha_pval < 0.05, na.rm = TRUE),
    n_pos_sig_05 = sum(alpha_pval < 0.05 & alpha > 0, na.rm = TRUE),
    n_neg_sig_05 = sum(alpha_pval < 0.05 & alpha < 0, na.rm = TRUE),
    avg_adj_r2 = mean(adj_r2, na.rm = TRUE)
  ),
  by = .(portfolio_set, variant, model)
][order(
  match(portfolio_set, chi_32_order),
  match(model, model_order)
)]

saveRDS(
  model_alpha_chi_32_vw_summary,
  "section10_model_alpha_chi_32_vw_summary.rds"
)

### 10.39E Requested SMQ Output ###############################################
smq_requested <- smq_summary_table[
  factor %in% c("SMQ_firm_VW", "SMQ_sic_VW", "SMQ_sic2_VW")
]

smq_requested_display <- copy(smq_requested)

smq_requested_display[, Mean_excess := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  100 * mean_excess, mean_stars, mean_tstat, mean_pval
)]

smq_requested_display[, CAPM_alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  100 * CAPM_alpha, CAPM_stars, CAPM_tstat, CAPM_pval
)]

smq_requested_display[, FF3_alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  100 * FF3_alpha, FF3_stars, FF3_tstat, FF3_pval
)]

smq_requested_display[, FF5_alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  100 * FF5_alpha, FF5_stars, FF5_tstat, FF5_pval
)]

smq_requested_display <- smq_requested_display[
  ,
  .(
    factor,
    Mean_excess,
    CAPM_alpha = CAPM_alpha_display,
    FF3_alpha  = FF3_alpha_display,
    FF5_alpha  = FF5_alpha_display
  )
][order(match(factor, c("SMQ_firm_VW", "SMQ_sic_VW", "SMQ_sic2_VW")))]

saveRDS(
  smq_requested_display,
  "section10_requested_smq_output_display.rds"
)

### 10.39F Requested Decile Portfolio Output ##################################
decile_requested_specs <- list(
  CHI_FIRM_DEC10_VW = list(dec_prefix = "CHI_FIRM_DD"),
  CHI_SIC_DEC10_VW  = list(dec_prefix = "CHI_SIC_DD"),
  CHI_SIC2_DEC10_VW = list(dec_prefix = "CHI_SIC2_DD")
)

build_requested_decile_detail <- function(panel, dec_prefix, lags = NW_LAGS) {
  dec_cols <- get_sorted_decile_cols(panel, dec_prefix)

  stopifnot(all(capm_rhs %in% names(panel)))
  stopifnot(all(ff3_rhs %in% names(panel)))
  stopifnot(all(ff5_rhs %in% names(panel)))
  stopifnot("RF" %in% names(panel))

  out <- rbindlist(lapply(seq_along(dec_cols), function(j) {
    col <- dec_cols[j]
    y <- panel[[col]] - panel[["RF"]]

    mean_res <- mean_nw(y, lags)
    capm_res <- run_ts_reg_nw(y, panel[, ..capm_rhs], lags)
    ff3_res <- run_ts_reg_nw(y, panel[, ..ff3_rhs], lags)
    ff5_res <- run_ts_reg_nw(y, panel[, ..ff5_rhs], lags)

    capm_a <- capm_res[term == "(Intercept)"][1]
    ff3_a <- ff3_res[term == "(Intercept)"][1]
    ff5_a <- ff5_res[term == "(Intercept)"][1]

    data.table(
      portfolio = paste0("P", j),
      mean_excess = mean_res$estimate,
      mean_tstat = mean_res$tstat,
      mean_pval = mean_res$pval,
      mean_stars = mean_res$stars,
      CAPM_alpha = capm_a$estimate,
      CAPM_tstat = capm_a$tstat,
      CAPM_pval = capm_a$pval,
      CAPM_stars = capm_a$stars,
      FF3_alpha = ff3_a$estimate,
      FF3_tstat = ff3_a$tstat,
      FF3_pval = ff3_a$pval,
      FF3_stars = ff3_a$stars,
      FF5_alpha = ff5_a$estimate,
      FF5_tstat = ff5_a$tstat,
      FF5_pval = ff5_a$pval,
      FF5_stars = ff5_a$stars
    )
  }))

  out
}

make_requested_decile_display <- function(detail_dt, family_name) {
  out <- copy(detail_dt)

  out[, Mean_excess := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * mean_excess, mean_stars, mean_tstat, mean_pval
  )]

  out[, CAPM := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * CAPM_alpha, CAPM_stars, CAPM_tstat, CAPM_pval
  )]

  out[, FF3 := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * FF3_alpha, FF3_stars, FF3_tstat, FF3_pval
  )]

  out[, FF5 := sprintf(
    "%.2f%s\n(t=%.2f, p=%.3f)",
    100 * FF5_alpha, FF5_stars, FF5_tstat, FF5_pval
  )]

  out[, family := family_name]
  out[, .(family, portfolio, Mean_excess, CAPM, FF3, FF5)]
}

make_requested_decile_sig_counts <- function(detail_dt, family_name) {
  data.table(
    family = family_name,
    n_mean_sig_05 = sum(detail_dt$mean_pval < 0.05, na.rm = TRUE),
    n_CAPM_sig_05 = sum(detail_dt$CAPM_pval < 0.05, na.rm = TRUE),
    n_FF3_sig_05 = sum(detail_dt$FF3_pval < 0.05, na.rm = TRUE),
    n_FF5_sig_05 = sum(detail_dt$FF5_pval < 0.05, na.rm = TRUE)
  )
}

requested_decile_detail_list <- list()
requested_decile_display_list <- list()
requested_decile_sigcount_list <- list()

for (nm in names(decile_requested_specs)) {
  detail_dt <- build_requested_decile_detail(
    panel = factor_test_panel,
    dec_prefix = decile_requested_specs[[nm]]$dec_prefix,
    lags = NW_LAGS
  )

  requested_decile_detail_list[[nm]] <- detail_dt
  requested_decile_display_list[[nm]] <- make_requested_decile_display(detail_dt, nm)
  requested_decile_sigcount_list[[nm]] <- make_requested_decile_sig_counts(detail_dt, nm)
}

requested_decile_detail_table <- rbindlist(
  lapply(names(requested_decile_detail_list), function(nm) {
    tmp <- copy(requested_decile_detail_list[[nm]])
    tmp[, family := nm]
    tmp
  }),
  fill = TRUE
)

requested_decile_display_table <- rbindlist(
  requested_decile_display_list,
  fill = TRUE
)

requested_decile_sigcount_table <- rbindlist(
  requested_decile_sigcount_list,
  fill = TRUE
)

saveRDS(
  requested_decile_detail_table,
  "section10_requested_decile_detail_table.rds"
)

saveRDS(
  requested_decile_display_table,
  "section10_requested_decile_display_table.rds"
)

saveRDS(
  requested_decile_sigcount_table,
  "section10_requested_decile_sigcount_table.rds"
)

### Identify significant decile alphas ########################################

requested_decile_detail_table <- readRDS("section10_requested_decile_detail_table.rds")
setDT(requested_decile_detail_table)

decile_alpha_long <- rbindlist(list(
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "CAPM",
      alpha = CAPM_alpha,
      tstat = CAPM_tstat,
      pval = CAPM_pval,
      stars = CAPM_stars
    )
  ],
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "FF3",
      alpha = FF3_alpha,
      tstat = FF3_tstat,
      pval = FF3_pval,
      stars = FF3_stars
    )
  ],
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "FF5",
      alpha = FF5_alpha,
      tstat = FF5_tstat,
      pval = FF5_pval,
      stars = FF5_stars
    )
  ]
), fill = TRUE)

decile_alpha_long[, construction := fcase(
  family == "CHI_FIRM_DEC10_VW", "Firm-level χ",
  family == "CHI_SIC_DEC10_VW", "Industry χ",
  family == "CHI_SIC2_DEC10_VW", "Major-group χ",
  default = family
)]

decile_alpha_long[, alpha_pct := 100 * alpha]

decile_alpha_long[, sig_level := fcase(
  pval < 0.01, "1%",
  pval < 0.05, "5%",
  pval < 0.10, "10%",
  default = "Not significant"
)]

sig_decile_alphas <- decile_alpha_long[
  pval < 0.10,
  .(
    construction,
    family,
    portfolio,
    model,
    alpha_pct,
    tstat,
    pval,
    stars,
    sig_level
  )
][order(construction, model, pval)]

print(round_dt(sig_decile_alphas, 4))

### Significant decile alphas by table threshold ##############################

make_sig_threshold <- function(threshold_label, cutoff) {
  out <- decile_alpha_long[
    pval < cutoff,
    .(
      threshold = threshold_label,
      construction,
      family,
      portfolio,
      model,
      alpha_pct,
      tstat,
      pval,
      stars
    )
  ]
  out[]
}

sig_decile_alphas_by_threshold <- rbindlist(list(
  make_sig_threshold("1% level", 0.01),
  make_sig_threshold("5% level", 0.05),
  make_sig_threshold("10% level", 0.10)
), fill = TRUE)

sig_decile_alphas_by_threshold <- sig_decile_alphas_by_threshold[
  order(threshold, construction, model, portfolio)
]

print(round_dt(sig_decile_alphas_by_threshold, 4))



### 10.39F1 Robustness Table: SMQ and Decile FF5 Alphas ########################
### Construction labels #######################################################
construction_order <- c(
  "Firm-level chi",
  "Industry chi",
  "Major-group chi"
)

construction_latex <- c(
  "Firm-level chi" = "Firm-level \\(\\chi\\)",
  "Industry chi" = "Industry \\(\\chi\\)",
  "Major-group chi" = "Major-group \\(\\chi\\)"
)

smq_factor_map <- data.table(
  factor = c("SMQ_firm_VW", "SMQ_sic_VW", "SMQ_sic2_VW"),
  construction = construction_order
)

decile_family_map <- data.table(
  family = c("CHI_FIRM_DEC10_VW", "CHI_SIC_DEC10_VW", "CHI_SIC2_DEC10_VW"),
  construction = construction_order
)

### Helper: LaTeX cell formatter ##############################################
latex_cell <- function(est, tstat, stars = "", digits_est = 2, digits_t = 2, scale = 100) {
  stars[is.na(stars)] <- ""

  out <- sprintf(
    paste0("\\makecell{%.", digits_est, "f%s \\\\ (%.", digits_t, "f)}"),
    scale * est,
    stars,
    tstat
  )

  out[is.na(est) | is.na(tstat)] <- "\\makecell{--}"
  out
}

latex_plain <- function(x, digits = 2) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), x))
}

### Panel A: SMQ factor properties ############################################
robustness_smq_panel <- merge(
  smq_summary_table[
    factor %in% smq_factor_map$factor
  ],
  smq_factor_map,
  by = "factor",
  all.x = TRUE
)

robustness_smq_panel[, construction := factor(
  construction,
  levels = construction_order
)]

robustness_smq_panel <- robustness_smq_panel[order(construction)]

robustness_smq_display_long <- rbindlist(list(
  robustness_smq_panel[
    ,
    .(
      row_label = "Mean excess return",
      construction,
      cell = latex_cell(mean_excess, mean_tstat, mean_stars)
    )
  ],
  robustness_smq_panel[
    ,
    .(
      row_label = "CAPM alpha",
      construction,
      cell = latex_cell(CAPM_alpha, CAPM_tstat, CAPM_stars)
    )
  ],
  robustness_smq_panel[
    ,
    .(
      row_label = "FF3 alpha",
      construction,
      cell = latex_cell(FF3_alpha, FF3_tstat, FF3_stars)
    )
  ],
  robustness_smq_panel[
    ,
    .(
      row_label = "FF5 alpha",
      construction,
      cell = latex_cell(FF5_alpha, FF5_tstat, FF5_stars)
    )
  ],
  robustness_smq_panel[
    ,
    .(
      row_label = "Monthly volatility",
      construction,
      cell = latex_plain(100 * vol_monthly, 2)
    )
  ],
  robustness_smq_panel[
    ,
    .(
      row_label = "Sharpe ratio",
      construction,
      cell = latex_plain(sharpe, 2)
    )
  ]
), fill = TRUE)

robustness_smq_display_long[, row_order := match(
  row_label,
  c(
    "Mean excess return",
    "CAPM alpha",
    "FF3 alpha",
    "FF5 alpha",
    "Monthly volatility",
    "Sharpe ratio"
  )
)]

robustness_smq_display <- dcast(
  robustness_smq_display_long,
  row_order + row_label ~ construction,
  value.var = "cell"
)[order(row_order)]

robustness_smq_display[, row_order := NULL]

### Panel B: FF5 decile alphas #################################################
robustness_ff5_decile_panel <- merge(
  requested_decile_detail_table[
    family %in% decile_family_map$family
  ],
  decile_family_map,
  by = "family",
  all.x = TRUE
)

robustness_ff5_decile_panel[, construction := factor(
  construction,
  levels = construction_order
)]

robustness_ff5_decile_panel[, decile_num := as.integer(gsub("[^0-9]", "", portfolio))]

robustness_ff5_decile_panel[, cell := latex_cell(
  FF5_alpha,
  FF5_tstat,
  FF5_stars
)]

robustness_ff5_decile_display <- dcast(
  robustness_ff5_decile_panel,
  decile_num + portfolio ~ construction,
  value.var = "cell"
)[order(decile_num)]

robustness_ff5_decile_display[, decile_num := NULL]

### Create LaTeX rows for direct copy-paste ###################################
make_latex_row <- function(row_name, firm_cell, industry_cell, major_cell) {
  paste0(
    row_name, "\n",
    "& ", firm_cell, "\n",
    "& ", industry_cell, "\n",
    "& ", major_cell, " \\\\"
  )
}

panel_a_rows <- robustness_smq_display[
  ,
  make_latex_row(
    row_label,
    get("Firm-level chi"),
    get("Industry chi"),
    get("Major-group chi")
  )
]

panel_b_rows <- robustness_ff5_decile_display[
  ,
  make_latex_row(
    portfolio,
    get("Firm-level chi"),
    get("Industry chi"),
    get("Major-group chi")
  )
]

robustness_table_latex_rows <- c(
  "\\multicolumn{4}{l}{\\textit{Panel A: SMQ factor properties}} \\\\",
  "\\midrule",
  panel_a_rows,
  "",
  "\\midrule",
  "\\multicolumn{4}{l}{\\textit{Panel B: FF5 alphas for \\(\\chi\\)-sorted decile portfolios}} \\\\",
  "\\midrule",
  panel_b_rows
)

### Save outputs ##############################################################
saveRDS(
  robustness_smq_display,
  "section10_robustness_smq_panel_display.rds"
)

saveRDS(
  robustness_ff5_decile_display,
  "section10_robustness_ff5_decile_alpha_display.rds"
)

saveRDS(
  robustness_table_latex_rows,
  "section10_robustness_table_latex_rows.rds"
)

writeLines(
  robustness_table_latex_rows,
  "section10_robustness_table_latex_rows.tex"
)

### Console output ############################################################
cat("\n------------------------------------------------------------\n")
cat("Robustness table Panel A: SMQ factor properties\n")
cat("------------------------------------------------------------\n")
print(robustness_smq_display)

cat("\n------------------------------------------------------------\n")
cat("Robustness table Panel B: FF5 decile alphas\n")
cat("------------------------------------------------------------\n")
print(robustness_ff5_decile_display)

cat("\n------------------------------------------------------------\n")
cat("LaTeX rows for robustness table\n")
cat("------------------------------------------------------------\n")
cat(paste(robustness_table_latex_rows, collapse = "\n"))
cat("\n")


### 10.39F2 FF5 Loadings of Alternative SMQ Specifications ####################
### Purpose:
### Regress each SMQ construction on the FF5 factors and extract the FF5 loadings.
### This shows how the different SMQ specifications relate to MKT, SMB, HML, RMW, and CMA.

smq_ff5_loading_specs <- data.table(
  factor = c("SMQ_firm_VW", "SMQ_sic_VW", "SMQ_sic2_VW"),
  construction = c("Firm-level chi", "Industry chi", "Major-group chi")
)

ff5_rhs <- c("MKT", "SMB", "HML", "RMW", "CMA")

make_smq_ff5_loading_table <- function(dt, smq_specs, rhs = ff5_rhs, lags = NW_LAGS) {
  out <- list()

  for (i in seq_len(nrow(smq_specs))) {
    smq_col <- smq_specs$factor[i]
    construction_name <- smq_specs$construction[i]

    stopifnot(smq_col %in% names(dt))
    stopifnot(all(rhs %in% names(dt)))

    reg <- run_ts_reg_nw_formula(
      as.formula(paste(smq_col, "~", paste(rhs, collapse = " + "))),
      dt,
      lags
    )

    if (is.null(reg)) next

    reg[, factor := smq_col]
    reg[, construction := construction_name]

    out[[smq_col]] <- reg
  }

  rbindlist(out, fill = TRUE)
}

smq_ff5_loading_table <- make_smq_ff5_loading_table(
  dt = panel,
  smq_specs = smq_ff5_loading_specs,
  rhs = ff5_rhs,
  lags = NW_LAGS
)

smq_ff5_loading_table[, term_label := fcase(
  term == "(Intercept)", "Alpha",
  term == "MKT", "MKT",
  term == "SMB", "SMB",
  term == "HML", "HML",
  term == "RMW", "RMW",
  term == "CMA", "CMA",
  default = term
)]

smq_ff5_loading_table[, term_order := match(
  term_label,
  c("Alpha", "MKT", "SMB", "HML", "RMW", "CMA")
)]

### Helper for LaTeX cells ####################################################

latex_cell_smq_loading <- function(est, tstat, stars = "", term_label, digits_est = 2, digits_t = 2) {
  stars[is.na(stars)] <- ""

  stars_latex <- fifelse(
    stars == "",
    "",
    paste0("$^{", stars, "}$")
  )

  scale_vec <- fifelse(term_label == "Alpha", 100, 1)

  out <- sprintf(
    paste0("\\makecell{%.", digits_est, "f%s \\\\ (%.", digits_t, "f)}"),
    scale_vec * est,
    stars_latex,
    tstat
  )

  out[is.na(est) | is.na(tstat)] <- "\\makecell{--}"
  out
}

smq_ff5_loading_display_long <- copy(smq_ff5_loading_table)

smq_ff5_loading_display_long[, cell := latex_cell_smq_loading(
  est = estimate,
  tstat = tstat,
  stars = stars,
  term_label = term_label
)]

smq_ff5_loading_display <- dcast(
  smq_ff5_loading_display_long,
  term_order + term_label ~ construction,
  value.var = "cell"
)[order(term_order)]

smq_ff5_loading_display[, term_order := NULL]

### Create LaTeX rows #########################################################

smq_ff5_loading_latex_rows <- smq_ff5_loading_display[
  ,
  make_latex_row(
    term_label,
    get("Firm-level chi"),
    get("Industry chi"),
    get("Major-group chi")
  )
]

smq_ff5_loading_latex_rows <- c(
  "\\multicolumn{4}{l}{\\textit{Panel C: FF5 loadings of alternative SMQ factors}} \\\\",
  "\\midrule",
  smq_ff5_loading_latex_rows
)

### Save outputs ##############################################################

saveRDS(
  smq_ff5_loading_table,
  "section10_smq_ff5_loading_table.rds"
)

saveRDS(
  smq_ff5_loading_display,
  "section10_smq_ff5_loading_display.rds"
)

saveRDS(
  smq_ff5_loading_latex_rows,
  "section10_smq_ff5_loading_latex_rows.rds"
)

writeLines(
  smq_ff5_loading_latex_rows,
  "section10_smq_ff5_loading_latex_rows.tex"
)

### Console output ############################################################

cat("\n------------------------------------------------------------\n")
cat("FF5 loadings of alternative SMQ specifications\n")
cat("------------------------------------------------------------\n")
print(smq_ff5_loading_display)

cat("\n------------------------------------------------------------\n")
cat("LaTeX rows for FF5 loadings of SMQ\n")
cat("------------------------------------------------------------\n")
cat(paste(smq_ff5_loading_latex_rows, collapse = "\n"))
cat("\n")





### 10.39G Requested HHI Summary Output #######################################
requested_hhi_count_summary <- factor_leg_hhi_summary_alt[
  Factor %in% c("SMQ_firm", "SMQ_sic", "SMQ_sic2")
][order(match(Factor, c("SMQ_firm", "SMQ_sic", "SMQ_sic2")), Leg)]

requested_hhi_weighted_summary <- factor_leg_hhi_weighted_summary_alt[
  Factor %in% c("SMQ_firm", "SMQ_sic", "SMQ_sic2")
][order(match(Factor, c("SMQ_firm", "SMQ_sic", "SMQ_sic2")), Leg)]

requested_hhi_overview <- merge(
  requested_hhi_count_summary,
  requested_hhi_weighted_summary[
    ,
    .(
      Factor,
      Leg,
      Avg_HHI_weighted = Avg_HHI,
      Median_HHI_weighted = Median_HHI,
      Min_HHI_weighted = Min_HHI,
      Max_HHI_weighted = Max_HHI,
      Avg_Effective_N_weighted = Avg_Effective_N
    )
  ],
  by = c("Factor", "Leg"),
  all = TRUE
)

saveRDS(
  requested_hhi_count_summary,
  "section10_requested_hhi_count_summary.rds"
)

saveRDS(
  requested_hhi_weighted_summary,
  "section10_requested_hhi_weighted_summary.rds"
)

saveRDS(
  requested_hhi_overview,
  "section10_requested_hhi_overview.rds"
)

### 10.39F1 Exact Significant Decile Portfolios ###############################

decile_alpha_long <- rbindlist(list(
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "Mean excess return",
      estimate = mean_excess,
      tstat = mean_tstat,
      pval = mean_pval,
      stars = mean_stars
    )
  ],
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "CAPM alpha",
      estimate = CAPM_alpha,
      tstat = CAPM_tstat,
      pval = CAPM_pval,
      stars = CAPM_stars
    )
  ],
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "FF3 alpha",
      estimate = FF3_alpha,
      tstat = FF3_tstat,
      pval = FF3_pval,
      stars = FF3_stars
    )
  ],
  requested_decile_detail_table[
    ,
    .(
      family,
      portfolio,
      model = "FF5 alpha",
      estimate = FF5_alpha,
      tstat = FF5_tstat,
      pval = FF5_pval,
      stars = FF5_stars
    )
  ]
), fill = TRUE)

decile_alpha_long[, construction := fcase(
  family == "CHI_FIRM_DEC10_VW", "Firm-level χ",
  family == "CHI_SIC_DEC10_VW", "Industry χ",
  family == "CHI_SIC2_DEC10_VW", "Major-group χ",
  default = family
)]

decile_alpha_long[, decile_num := as.integer(gsub("[^0-9]", "", portfolio))]
decile_alpha_long[, estimate_pct := 100 * estimate]

decile_alpha_long[, sig_level := fcase(
  pval < 0.01, "1%",
  pval < 0.05, "5%",
  pval < 0.10, "10%",
  default = "Not significant"
)]

decile_alpha_long[, sign := fcase(
  estimate > 0, "Positive",
  estimate < 0, "Negative",
  default = "Zero"
)]

decile_alpha_long[, chi_region := fcase(
  decile_num <= 3, "Low χ",
  decile_num >= 8, "High χ",
  default = "Middle χ"
)]

significant_decile_portfolios <- decile_alpha_long[
  pval < 0.10,
  .(
    construction,
    family,
    portfolio,
    decile_num,
    chi_region,
    model,
    estimate_pct,
    tstat,
    pval,
    stars,
    sig_level,
    sign
  )
][order(construction, model, decile_num)]

print(round_dt(significant_decile_portfolios, 4))

saveRDS(
  decile_alpha_long,
  "section10_decile_alpha_long.rds"
)

saveRDS(
  significant_decile_portfolios,
  "section10_significant_decile_portfolios.rds"
)

### 10.39H Consolidated Master Output #########################################
section10_requested_capm_and_compact_outputs <- list(
  model_grs_chi_25_vw = model_grs_chi_25_vw,
  model_grs_chi_25_vw_display = model_grs_chi_25_vw_display,
  model_alpha_chi_25_vw = model_alpha_chi_25_vw,
  model_alpha_chi_25_vw_display = model_alpha_chi_25_vw_display,
  model_alpha_chi_25_vw_summary = model_alpha_chi_25_vw_summary,
  model_grs_chi_32_vw = model_grs_chi_32_vw,
  model_grs_chi_32_vw_display = model_grs_chi_32_vw_display,
  model_alpha_chi_32_vw = model_alpha_chi_32_vw,
  model_alpha_chi_32_vw_display = model_alpha_chi_32_vw_display,
  model_alpha_chi_32_vw_summary = model_alpha_chi_32_vw_summary,
  smq_requested_display = smq_requested_display,
  requested_decile_detail_table = requested_decile_detail_table,
  requested_decile_display_table = requested_decile_display_table,
  requested_decile_sigcount_table = requested_decile_sigcount_table,
  requested_hhi_count_summary = requested_hhi_count_summary,
  requested_hhi_weighted_summary = requested_hhi_weighted_summary,
  requested_hhi_overview = requested_hhi_overview
)

saveRDS(
  section10_requested_capm_and_compact_outputs,
  "section10_requested_capm_and_compact_outputs.rds"
)

### 10.39I Total Console Summary ##############################################
## SMQ Factor Test
cat("------------------------------------------------------------\n")
cat("SMQ summary (firm VW, SIC VW, SIC2 VW)\n")
cat("------------------------------------------------------------\n")
print(smq_requested_display)

cat("------------------------------------------------------------\n")
cat("Full SMQ summary table (all underlying statistics)\n")
cat("------------------------------------------------------------\n")
print(round_dt(smq_summary_table, 4))

## Chi-decile tests
cat("\n------------------------------------------------------------\n")
cat("Decile output: significance counts at 5%\n")
cat("------------------------------------------------------------\n")
print(round_dt(requested_decile_sigcount_table, 4))

## Chi decile tables
cat("\n------------------------------------------------------------\n")
cat("Requested decile detail table\n")
cat("------------------------------------------------------------\n")
print(requested_decile_display_table)

## Industry concentration
cat("\n------------------------------------------------------------\n")
cat("HHI overview summary\n")
cat("------------------------------------------------------------\n")
print(round_dt(requested_hhi_overview, 4))

## Test Portfolios 5x5 Sorts
# GRS Test
cat("\n------------------------------------------------------------\n")
cat("GRS results: 25-portfolio chi blocks (VW only)\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_grs_chi_25_vw_display, 4))

## Number of significant alphas
cat("\n------------------------------------------------------------\n")
cat("Alpha summary: 25-portfolio chi blocks (VW only)\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_alpha_chi_25_vw_summary, 4))

# Individual models granular results
cat("\n------------------------------------------------------------\n")
cat("Individual alphas: 25-portfolio chi blocks\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_FIRM_25" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_FIRM_25" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_FIRM_25" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_FIRM_25" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC_25" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC_25" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC_25" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC_25" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC2_25" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC2_25" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC2_25" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_25_vw_display[portfolio_set == "CHI_SIC2_25" & model == "FF5_SMQ"], 4))

## Test Portfolios 2x4x4 Sorts
# GRS Test
cat("\n------------------------------------------------------------\n")
cat("GRS results: 32-portfolio chi families (VW only)\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_grs_chi_32_vw_display, 4))

## Number of significant alphas
cat("\n------------------------------------------------------------\n")
cat("Alpha summary: 32-portfolio chi families (VW only)\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_alpha_chi_32_vw_summary, 4))

# Indidvidual models granular results
cat("\n------------------------------------------------------------\n")
cat("Individual alphas: 32-portfolio chi families\n")
cat("------------------------------------------------------------\n")
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_FIRM_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_FIRM_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_FIRM_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_FIRM_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_FIRM_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_FIRM_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_FIRM_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_FIRM_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_FIRM_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_FIRM_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_FIRM_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_FIRM_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC2_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC2_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC2_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "OP_CHI_SIC2_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC2_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC2_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC2_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "BM_CHI_SIC2_32" & model == "FF5_SMQ"], 4))

print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC2_32" & model == "CAPM"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC2_32" & model == "FF3"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC2_32" & model == "FF5"], 4))
print(round_dt(model_alpha_chi_32_vw_display[portfolio_set == "INV_CHI_SIC2_32" & model == "FF5_SMQ"], 4))

## Simple count of significant alphas on decile portfolios
smq_panel_counts <- requested_decile_detail_table[
  ,
  .(
    mean_1  = sum(mean_pval < 0.01, na.rm = TRUE),
    capm_1  = sum(CAPM_pval < 0.01, na.rm = TRUE),
    ff3_1   = sum(FF3_pval < 0.01, na.rm = TRUE),
    ff5_1   = sum(FF5_pval < 0.01, na.rm = TRUE),
    mean_5  = sum(mean_pval < 0.05, na.rm = TRUE),
    capm_5  = sum(CAPM_pval < 0.05, na.rm = TRUE),
    ff3_5   = sum(FF3_pval < 0.05, na.rm = TRUE),
    ff5_5   = sum(FF5_pval < 0.05, na.rm = TRUE),
    mean_10 = sum(mean_pval < 0.10, na.rm = TRUE),
    capm_10 = sum(CAPM_pval < 0.10, na.rm = TRUE),
    ff3_10  = sum(FF3_pval < 0.10, na.rm = TRUE),
    ff5_10  = sum(FF5_pval < 0.10, na.rm = TRUE)
  ),
  by = .(family)
]

smq_panel_counts[, family := factor(
  family,
  levels = c("CHI_FIRM_DEC10_VW", "CHI_SIC_DEC10_VW", "CHI_SIC2_DEC10_VW"),
  labels = c("Firm-level χ", "Industry χ", "Major-group χ")
)]

setnames(smq_panel_counts, "family", "construction")

print(smq_panel_counts)

## Simple count of significant alphas on 5x5 test portfolios and 2x4x4 test portfolios
count_sig_levels <- function(dt) {
  data.table(
    sig_1  = sum(dt$alpha_pval < 0.01, na.rm = TRUE),
    sig_5  = sum(dt$alpha_pval < 0.05, na.rm = TRUE),
    sig_10 = sum(dt$alpha_pval < 0.10, na.rm = TRUE)
  )
}

chi_alpha_counts_all <- rbindlist(list(
  model_alpha_chi_25_vw_display,
  model_alpha_chi_32_vw_display
), fill = TRUE)

chi_alpha_counts_summary <- chi_alpha_counts_all[
  model %in% c("CAPM", "FF3", "FF5", "FF5_SMQ"),
  .(
    sig_1  = sum(alpha_pval < 0.01, na.rm = TRUE),
    sig_5  = sum(alpha_pval < 0.05, na.rm = TRUE),
    sig_10 = sum(alpha_pval < 0.10, na.rm = TRUE)
  ),
  by = .(portfolio_set, model)
][order(portfolio_set, factor(model, levels = c("CAPM", "FF3", "FF5", "FF5_SMQ")))]

print(chi_alpha_counts_summary)



### 10.40 Overlap between INV_CHI S44 and OP_INV S14 ###########################
### 1. Inspect available portfolio-label columns ##############################
port_cols <- grep("port|OP_INV|INV_CHI|CHI_INV", names(labelled_data), value = TRUE)
print(port_cols)

### 2. Choose relevant portfolio-label columns ###############################
## OP_INV column: small OP=1, Inv=4 is S14
op_inv_col <- "port_OP_INV_2x4x4"

## INV_CHI column: small Inv=4, Chi=4 is S44
## Change this if you want firm-level or major-group chi instead.
## Likely candidates may be:
## "port_INV_CHI_FIRM_2x4x4"
## "port_INV_CHI_SIC_2x4x4"
## "port_INV_CHI_SIC2_2x4x4"
## "port_INV_CHI_2x4x4"

inv_chi_candidates <- c(
  "port_INV_CHI_SIC_2x4x4",
  "port_INV_CHI_2x4x4",
  "port_INV_CHI_FIRM_2x4x4",
  "port_INV_CHI_SIC2_2x4x4",
  "port_INV_CHI_SIC",
  "port_INV_CHI_FIRM",
  "port_INV_CHI_SIC2"
)

inv_chi_col <- inv_chi_candidates[inv_chi_candidates %in% names(labelled_data)][1]

if (is.na(inv_chi_col)) {
  stop("No INV_CHI portfolio-label column found. Inspect `port_cols` above and set `inv_chi_col` manually.")
}

if (!op_inv_col %in% names(labelled_data)) {
  stop("OP_INV portfolio-label column not found. Inspect `port_cols` above and set `op_inv_col` manually.")
}

cat("Using OP_INV column: ", op_inv_col, "\n")
cat("Using INV_CHI column:", inv_chi_col, "\n")


### 3. Helper function: robust label match ####################################
## Works if labels are stored as "S44" or as something ending in "_S44".
match_port_label <- function(x, target_label) {
  x <- as.character(x)
  !is.na(x) & (
    x == target_label |
      grepl(paste0("(^|_)", target_label, "$"), x)
  )
}


### 4. Main overlap function ##################################################
make_overlap_diagnostics <- function(
    dt,
    port_a_col,
    port_a_label,
    port_b_col,
    port_b_label,
    id_col = "permno",
    weight_candidates = c("me_ff_w", "me_clean", "me_ff", "me_june")) {
  stopifnot(port_a_col %in% names(dt))
  stopifnot(port_b_col %in% names(dt))
  stopifnot(id_col %in% names(dt))
  stopifnot("mdate" %in% names(dt))

  x <- copy(dt)
  setDT(x)
  x[, mdate := as.Date(mdate)]

  weight_col <- weight_candidates[weight_candidates %in% names(x)][1]
  if (is.na(weight_col)) {
    weight_col <- NULL
    x[, overlap_weight := 1]
  } else {
    x[, overlap_weight := get(weight_col)]
    x[is.na(overlap_weight) | !is.finite(overlap_weight) | overlap_weight <= 0, overlap_weight := NA_real_]
  }

  x[, in_A := match_port_label(get(port_a_col), port_a_label)]
  x[, in_B := match_port_label(get(port_b_col), port_b_label)]
  x[, in_both := in_A & in_B]
  x[, in_union := in_A | in_B]

  ### Monthly overlap
  monthly <- x[
    in_union == TRUE,
    .(
      N_A     = uniqueN(get(id_col)[in_A]),
      N_B     = uniqueN(get(id_col)[in_B]),
      N_both  = uniqueN(get(id_col)[in_both]),
      N_union = uniqueN(get(id_col)[in_union]),
      W_A     = sum(overlap_weight[in_A], na.rm = TRUE),
      W_B     = sum(overlap_weight[in_B], na.rm = TRUE),
      W_both  = sum(overlap_weight[in_both], na.rm = TRUE),
      W_union = sum(overlap_weight[in_union], na.rm = TRUE)
    ),
    by = mdate
  ]

  monthly[, overlap_A_share := N_both / N_A]
  monthly[, overlap_B_share := N_both / N_B]
  monthly[, jaccard_share := N_both / N_union]

  monthly[, overlap_A_weight_share := W_both / W_A]
  monthly[, overlap_B_weight_share := W_both / W_B]
  monthly[, jaccard_weight_share := W_both / W_union]

  ### Formation-year overlap, if ffyear exists
  if ("ffyear" %in% names(x)) {
    form <- unique(
      x[
        in_union == TRUE,
        .(
          ffyear,
          id = get(id_col),
          in_A,
          in_B,
          in_both,
          in_union,
          overlap_weight
        )
      ],
      by = c("ffyear", "id", "in_A", "in_B")
    )

    formation <- form[
      ,
      .(
        N_A     = uniqueN(id[in_A]),
        N_B     = uniqueN(id[in_B]),
        N_both  = uniqueN(id[in_both]),
        N_union = uniqueN(id[in_union]),
        W_A     = sum(overlap_weight[in_A], na.rm = TRUE),
        W_B     = sum(overlap_weight[in_B], na.rm = TRUE),
        W_both  = sum(overlap_weight[in_both], na.rm = TRUE),
        W_union = sum(overlap_weight[in_union], na.rm = TRUE)
      ),
      by = ffyear
    ]

    formation[, overlap_A_share := N_both / N_A]
    formation[, overlap_B_share := N_both / N_B]
    formation[, jaccard_share := N_both / N_union]

    formation[, overlap_A_weight_share := W_both / W_A]
    formation[, overlap_B_weight_share := W_both / W_B]
    formation[, jaccard_weight_share := W_both / W_union]
  } else {
    formation <- NULL
  }

  ### Full-sample summary
  summary_monthly <- monthly[
    ,
    .(
      avg_N_A = mean(N_A, na.rm = TRUE),
      avg_N_B = mean(N_B, na.rm = TRUE),
      avg_N_both = mean(N_both, na.rm = TRUE),
      avg_N_union = mean(N_union, na.rm = TRUE),
      avg_overlap_A_share = mean(overlap_A_share, na.rm = TRUE),
      avg_overlap_B_share = mean(overlap_B_share, na.rm = TRUE),
      avg_jaccard_share = mean(jaccard_share, na.rm = TRUE),
      avg_overlap_A_weight_share = mean(overlap_A_weight_share, na.rm = TRUE),
      avg_overlap_B_weight_share = mean(overlap_B_weight_share, na.rm = TRUE),
      avg_jaccard_weight_share = mean(jaccard_weight_share, na.rm = TRUE)
    )
  ]

  list(
    weight_col_used = weight_col,
    monthly = monthly,
    formation = formation,
    summary_monthly = summary_monthly
  )
}


### 5. Run overlap test #######################################################
overlap_INVCHI_S44_OPINV_S14 <- make_overlap_diagnostics(
  dt = labelled_data[
    mdate >= START_DATE &
      mdate <= END_DATE
  ],
  port_a_col = inv_chi_col,
  port_a_label = "S44", # Small, Inv = 4, Chi = 4
  port_b_col = op_inv_col,
  port_b_label = "S14" # Small, OP = 1, Inv = 4
)

cat("\nWeight column used:\n")
print(overlap_INVCHI_S44_OPINV_S14$weight_col_used)

cat("\nMonthly overlap summary:\n")
print(round_dt(overlap_INVCHI_S44_OPINV_S14$summary_monthly, 4))

cat("\nFirst rows of monthly overlap:\n")
print(round_dt(head(overlap_INVCHI_S44_OPINV_S14$monthly, 10), 4))

if (!is.null(overlap_INVCHI_S44_OPINV_S14$formation)) {
  cat("\nFirst rows of formation-year overlap:\n")
  print(round_dt(head(overlap_INVCHI_S44_OPINV_S14$formation, 10), 4))
}

saveRDS(
  overlap_INVCHI_S44_OPINV_S14,
  "section10_overlap_INVCHI_S44_OPINV_S14.rds"
)




### 10.41 Robustness of Small-Firm Chi Loading Channel #########################
### 10.41A Helper functions ####################################################

get_portfolio_series_2x4x4 <- function(panel, prefix, size_flag, char_bucket, chi_bucket) {
  col <- paste0(prefix, size_flag, char_bucket, chi_bucket)

  if (!col %in% names(panel)) {
    stop("Portfolio column not found: ", col)
  }

  panel[[col]]
}

run_spread_loading_reg <- function(dt, y, rhs, meta, lags = NW_LAGS) {
  reg <- run_ts_reg_nw(
    y = y,
    X = dt[, ..rhs],
    lags = lags
  )

  if (is.null(reg)) {
    return(NULL)
  }

  out <- copy(reg)

  for (nm in names(meta)) {
    out[, (nm) := meta[[nm]]]
  }

  out[]
}

make_chi_loading_spread_tests_one_family <- function(
    dt,
    family_name,
    prefix,
    rhs = c("MKT", "SMB", "HML", "RMW", "CMA"),
    period_name = "Full sample",
    lags = NW_LAGS) {
  out <- list()

  for (size_flag in c("S", "B")) {
    spread_mat <- matrix(NA_real_, nrow = nrow(dt), ncol = 4)

    for (ii in 1:4) {
      low_series <- get_portfolio_series_2x4x4(
        panel = dt,
        prefix = prefix,
        size_flag = size_flag,
        char_bucket = ii,
        chi_bucket = 1
      )

      high_series <- get_portfolio_series_2x4x4(
        panel = dt,
        prefix = prefix,
        size_flag = size_flag,
        char_bucket = ii,
        chi_bucket = 4
      )

      y_spread <- high_series - low_series
      spread_mat[, ii] <- y_spread

      reg <- run_spread_loading_reg(
        dt = dt,
        y = y_spread,
        rhs = rhs,
        meta = list(
          period = period_name,
          family = family_name,
          spread_type = "Within size: high chi minus low chi",
          size_bucket = size_flag,
          char_bucket = paste0("Q", ii)
        ),
        lags = lags
      )

      out[[paste(family_name, period_name, size_flag, ii, "HL", sep = "_")]] <- reg
    }

    y_avg <- rowMeans(spread_mat, na.rm = FALSE)

    reg_avg <- run_spread_loading_reg(
      dt = dt,
      y = y_avg,
      rhs = rhs,
      meta = list(
        period = period_name,
        family = family_name,
        spread_type = "Within size: average high chi minus low chi",
        size_bucket = size_flag,
        char_bucket = "Average"
      ),
      lags = lags
    )

    out[[paste(family_name, period_name, size_flag, "Average_HL", sep = "_")]] <- reg_avg
  }

  for (ii in 1:4) {
    s_low <- get_portfolio_series_2x4x4(dt, prefix, "S", ii, 1)
    s_high <- get_portfolio_series_2x4x4(dt, prefix, "S", ii, 4)
    b_low <- get_portfolio_series_2x4x4(dt, prefix, "B", ii, 1)
    b_high <- get_portfolio_series_2x4x4(dt, prefix, "B", ii, 4)

    y_double <- (s_high - s_low) - (b_high - b_low)

    reg_double <- run_spread_loading_reg(
      dt = dt,
      y = y_double,
      rhs = rhs,
      meta = list(
        period = period_name,
        family = family_name,
        spread_type = "Double difference: small chi spread minus big chi spread",
        size_bucket = "S_minus_B",
        char_bucket = paste0("Q", ii)
      ),
      lags = lags
    )

    out[[paste(family_name, period_name, "Double", ii, sep = "_")]] <- reg_double
  }

  s_spreads <- matrix(NA_real_, nrow = nrow(dt), ncol = 4)
  b_spreads <- matrix(NA_real_, nrow = nrow(dt), ncol = 4)

  for (ii in 1:4) {
    s_spreads[, ii] <- get_portfolio_series_2x4x4(dt, prefix, "S", ii, 4) -
      get_portfolio_series_2x4x4(dt, prefix, "S", ii, 1)

    b_spreads[, ii] <- get_portfolio_series_2x4x4(dt, prefix, "B", ii, 4) -
      get_portfolio_series_2x4x4(dt, prefix, "B", ii, 1)
  }

  y_double_avg <- rowMeans(s_spreads, na.rm = FALSE) - rowMeans(b_spreads, na.rm = FALSE)

  reg_double_avg <- run_spread_loading_reg(
    dt = dt,
    y = y_double_avg,
    rhs = rhs,
    meta = list(
      period = period_name,
      family = family_name,
      spread_type = "Double difference: average small chi spread minus average big chi spread",
      size_bucket = "S_minus_B",
      char_bucket = "Average"
    ),
    lags = lags
  )

  out[[paste(family_name, period_name, "Double_Average", sep = "_")]] <- reg_double_avg

  rbindlist(out, fill = TRUE)
}

make_chi_loading_spread_tests_all <- function(
    panel,
    family_specs,
    periods,
    rhs = c("MKT", "SMB", "HML", "RMW", "CMA"),
    lags = NW_LAGS) {
  out <- list()

  for (period_name in names(periods)) {
    period_spec <- periods[[period_name]]

    dt_p <- panel[
      mdate >= period_spec$start &
        mdate <= period_spec$end
    ]

    for (family_name in names(family_specs)) {
      prefix <- family_specs[[family_name]]$prefix

      reg_dt <- make_chi_loading_spread_tests_one_family(
        dt = dt_p,
        family_name = family_name,
        prefix = prefix,
        rhs = rhs,
        period_name = period_name,
        lags = lags
      )

      out[[paste(period_name, family_name, sep = "_")]] <- reg_dt
    }
  }

  rbindlist(out, fill = TRUE)
}

### 10.41B Run loading-spread robustness tests ################################

chi_loading_periods <- list(
  Full = list(start = START_DATE, end = END_DATE),
  Early = list(start = as.Date("1982-07-01"), end = as.Date("2002-12-01")),
  Late = list(start = as.Date("2003-01-01"), end = END_DATE)
)

chi_loading_family_specs <- chi_32_specs

chi_loading_spread_tests <- make_chi_loading_spread_tests_all(
  panel = factor_test_panel,
  family_specs = chi_loading_family_specs,
  periods = chi_loading_periods,
  rhs = c("MKT", "SMB", "HML", "RMW", "CMA"),
  lags = NW_LAGS
)

chi_loading_spread_tests[, sort_family := fcase(
  grepl("^OP_", family), "Size-OP-chi",
  grepl("^BM_", family), "Size-BM-chi",
  grepl("^INV_", family), "Size-Inv-chi",
  default = family
)]

chi_loading_spread_tests[, chi_construction := fcase(
  grepl("FIRM", family), "Firm-level chi",
  grepl("SIC2", family), "Major-group chi",
  grepl("SIC", family), "Industry chi",
  default = family
)]

chi_loading_spread_tests[, estimate_pct := 100 * estimate]

### 10.41C Core HML and RMW loading-spread table ##############################

chi_loading_core <- chi_loading_spread_tests[
  term %in% c("(Intercept)", "HML", "RMW") &
    spread_type %in% c(
      "Within size: high chi minus low chi",
      "Within size: average high chi minus low chi",
      "Double difference: small chi spread minus big chi spread",
      "Double difference: average small chi spread minus average big chi spread"
    ),
  .(
    period,
    family,
    sort_family,
    chi_construction,
    spread_type,
    size_bucket,
    char_bucket,
    term,
    estimate,
    estimate_pct,
    tstat,
    pval,
    stars,
    adj_r2
  )
][order(period, chi_construction, sort_family, spread_type, size_bucket, char_bucket, term)]

saveRDS(
  chi_loading_spread_tests,
  "section10_chi_loading_spread_tests_full.rds"
)

saveRDS(
  chi_loading_core,
  "section10_chi_loading_core_hml_rmw.rds"
)

### 10.41D Sign-based robustness summary ######################################

chi_loading_sign_summary <- chi_loading_core[
  term %in% c("HML", "RMW") &
    spread_type == "Within size: high chi minus low chi",
  .(
    n_cells = .N,
    mean_loading_spread = mean(estimate, na.rm = TRUE),
    median_loading_spread = median(estimate, na.rm = TRUE),
    share_negative = mean(estimate < 0, na.rm = TRUE),
    n_negative = sum(estimate < 0, na.rm = TRUE),
    n_negative_sig_10 = sum(estimate < 0 & pval < 0.10, na.rm = TRUE),
    n_negative_sig_05 = sum(estimate < 0 & pval < 0.05, na.rm = TRUE)
  ),
  by = .(period, sort_family, chi_construction, size_bucket, term)
][order(period, sort_family, chi_construction, size_bucket, term)]

saveRDS(
  chi_loading_sign_summary,
  "section10_chi_loading_sign_summary.rds"
)

### 10.41E Small-minus-big robustness summary #################################

chi_loading_double_diff_summary <- chi_loading_core[
  term %in% c("HML", "RMW") &
    spread_type == "Double difference: small chi spread minus big chi spread",
  .(
    n_cells = .N,
    mean_double_diff = mean(estimate, na.rm = TRUE),
    median_double_diff = median(estimate, na.rm = TRUE),
    share_negative = mean(estimate < 0, na.rm = TRUE),
    n_negative = sum(estimate < 0, na.rm = TRUE),
    n_negative_sig_10 = sum(estimate < 0 & pval < 0.10, na.rm = TRUE),
    n_negative_sig_05 = sum(estimate < 0 & pval < 0.05, na.rm = TRUE)
  ),
  by = .(period, sort_family, chi_construction, term)
][order(period, sort_family, chi_construction, term)]

saveRDS(
  chi_loading_double_diff_summary,
  "section10_chi_loading_double_diff_summary.rds"
)

### 10.41F Average-spread display table #######################################

chi_loading_average_display <- chi_loading_core[
  char_bucket == "Average" &
    term %in% c("(Intercept)", "HML", "RMW"),
  .(
    period,
    family,
    sort_family,
    chi_construction,
    spread_type,
    size_bucket,
    term,
    estimate,
    estimate_pct,
    tstat,
    pval,
    stars,
    adj_r2
  )
][order(period, chi_construction, sort_family, spread_type, size_bucket, term)]

saveRDS(
  chi_loading_average_display,
  "section10_chi_loading_average_display.rds"
)

### 10.41G Console output #####################################################

cat("\n============================================================\n")
cat("ROBUSTNESS: SMALL-FIRM CHI LOADING CHANNEL\n")
cat("============================================================\n\n")

cat("Interpretation guide:\n")
cat("- Negative HML high-minus-low chi spread: high-chi portfolios are more growth-like.\n")
cat("- Negative RMW high-minus-low chi spread: high-chi portfolios are more weak-profitability-like.\n")
cat("- More negative RMW for small firms supports the small-firm profitability channel.\n")
cat("- Negative small-minus-big double difference supports a stronger small-firm chi-loading effect.\n\n")

cat("Average high-minus-low chi spread regressions:\n")
print(round_dt(chi_loading_average_display, 4))

cat("\nSign-based summary across individual 4x4 cells:\n")
print(round_dt(chi_loading_sign_summary, 4))

cat("\nSmall-minus-big double-difference summary:\n")
print(round_dt(chi_loading_double_diff_summary, 4))



### 10.42 FF5 Alphas of High-Minus-Low Chi Spread Portfolios ##################
### 10.42A Helper functions ###################################################
format_alpha_latex_cell <- function(alpha, tstat, stars, digits_est = 2, digits_t = 2, scale = 100) {
  if (is.na(alpha) | is.na(tstat)) {
    return("")
  }

  star_latex <- ifelse(
    is.na(stars) | stars == "",
    "",
    paste0("^{", stars, "}")
  )

  sprintf(
    paste0("\\makecell{$%.", digits_est, "f%s$ \\\\ $(%.", digits_t, "f)$}"),
    scale * alpha,
    star_latex,
    tstat
  )
}

run_ff5_alpha_spread <- function(dt, y, meta, lags = NW_LAGS) {
  reg <- run_ts_reg_nw(
    y = y,
    X = dt[, .(MKT, SMB, HML, RMW, CMA)],
    lags = lags
  )

  if (is.null(reg)) {
    out <- data.table(
      alpha = NA_real_,
      alpha_tstat = NA_real_,
      alpha_pval = NA_real_,
      alpha_stars = "",
      adj_r2 = NA_real_
    )
  } else {
    a <- reg[term == "(Intercept)"][1]

    out <- data.table(
      alpha = a$estimate,
      alpha_tstat = a$tstat,
      alpha_pval = a$pval,
      alpha_stars = a$stars,
      adj_r2 = a$adj_r2
    )
  }

  for (nm in names(meta)) {
    out[, (nm) := meta[[nm]]]
  }

  out[]
}

get_portfolio_return <- function(dt, col) {
  if (!col %in% names(dt)) {
    stop("Portfolio column not found: ", col)
  }

  dt[[col]]
}

### 10.42B Size-chi 25 portfolios: C5 minus C1 ###############################

make_size_chi_hl_ff5_alphas <- function(
    dt,
    construction_name,
    prefix,
    portfolio_set,
    lags = NW_LAGS) {
  out <- list()
  size_spreads <- list()

  for (ss in 1:5) {
    low_col <- paste0(prefix, "S", ss, "C1")
    high_col <- paste0(prefix, "S", ss, "C5")

    y <- get_portfolio_return(dt, high_col) - get_portfolio_return(dt, low_col)
    size_spreads[[paste0("S", ss)]] <- y

    out[[paste0("S", ss)]] <- run_ff5_alpha_spread(
      dt = dt,
      y = y,
      meta = list(
        portfolio_set = portfolio_set,
        construction = construction_name,
        spread_family = "25 portfolios: Size-chi",
        spread_row = paste0("S", ss),
        spread_definition = paste0("S", ss, ": C5 minus C1"),
        compact_row = NA_character_
      ),
      lags = lags
    )
  }

  y_small <- size_spreads[["S1"]]
  y_big <- size_spreads[["S5"]]
  y_small_big <- y_small - y_big

  out[["Small"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_small,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = "25 portfolios: Size-chi",
      spread_row = "Small",
      spread_definition = "Small: S1 C5 minus S1 C1",
      compact_row = "Small"
    ),
    lags = lags
  )

  out[["Big"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_big,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = "25 portfolios: Size-chi",
      spread_row = "Big",
      spread_definition = "Big: S5 C5 minus S5 C1",
      compact_row = "Big"
    ),
    lags = lags
  )

  out[["Small_Big"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_small_big,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = "25 portfolios: Size-chi",
      spread_row = "Small-Big",
      spread_definition = "(S1 C5 minus S1 C1) minus (S5 C5 minus S5 C1)",
      compact_row = "Small-Big"
    ),
    lags = lags
  )

  rbindlist(out, fill = TRUE)
}

size_chi_hl_specs <- list(
  FIRM = list(
    construction_name = "Firm",
    prefix = "CHI_FIRM_",
    portfolio_set = "Size_chi"
  ),
  INDUSTRY = list(
    construction_name = "Industry",
    prefix = "CHI_SIC_",
    portfolio_set = "Size_chi"
  ),
  MAJOR_GROUP = list(
    construction_name = "Major group",
    prefix = "CHI_SIC2_",
    portfolio_set = "Size_chi"
  )
)

size_chi_hl_ff5_alpha <- rbindlist(
  lapply(size_chi_hl_specs, function(sp) {
    make_size_chi_hl_ff5_alphas(
      dt = factor_test_panel,
      construction_name = sp$construction_name,
      prefix = sp$prefix,
      portfolio_set = sp$portfolio_set,
      lags = NW_LAGS
    )
  }),
  fill = TRUE
)

### 10.42C 32 portfolios: chi Q4 minus chi Q1 averaged across characteristic buckets

make_32_chi_hl_ff5_alphas <- function(
    dt,
    construction_name,
    prefix,
    portfolio_set,
    spread_family,
    lags = NW_LAGS) {
  out <- list()

  spread_mat_s <- matrix(NA_real_, nrow = nrow(dt), ncol = 4)
  spread_mat_b <- matrix(NA_real_, nrow = nrow(dt), ncol = 4)

  for (ii in 1:4) {
    low_s <- paste0(prefix, "S", ii, "1")
    high_s <- paste0(prefix, "S", ii, "4")
    low_b <- paste0(prefix, "B", ii, "1")
    high_b <- paste0(prefix, "B", ii, "4")

    y_s <- get_portfolio_return(dt, high_s) - get_portfolio_return(dt, low_s)
    y_b <- get_portfolio_return(dt, high_b) - get_portfolio_return(dt, low_b)

    spread_mat_s[, ii] <- y_s
    spread_mat_b[, ii] <- y_b

    out[[paste0("S_Q", ii)]] <- run_ff5_alpha_spread(
      dt = dt,
      y = y_s,
      meta = list(
        portfolio_set = portfolio_set,
        construction = construction_name,
        spread_family = spread_family,
        spread_row = paste0("Small Q", ii),
        spread_definition = paste0("Small, characteristic Q", ii, ": chi Q4 minus chi Q1"),
        compact_row = NA_character_
      ),
      lags = lags
    )

    out[[paste0("B_Q", ii)]] <- run_ff5_alpha_spread(
      dt = dt,
      y = y_b,
      meta = list(
        portfolio_set = portfolio_set,
        construction = construction_name,
        spread_family = spread_family,
        spread_row = paste0("Big Q", ii),
        spread_definition = paste0("Big, characteristic Q", ii, ": chi Q4 minus chi Q1"),
        compact_row = NA_character_
      ),
      lags = lags
    )
  }

  y_small <- rowMeans(spread_mat_s, na.rm = FALSE)
  y_big <- rowMeans(spread_mat_b, na.rm = FALSE)
  y_small_big <- y_small - y_big

  out[["Small_average"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_small,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = spread_family,
      spread_row = "Small",
      spread_definition = "Small average: average chi Q4 minus chi Q1 across characteristic quartiles",
      compact_row = "Small"
    ),
    lags = lags
  )

  out[["Big_average"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_big,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = spread_family,
      spread_row = "Big",
      spread_definition = "Big average: average chi Q4 minus chi Q1 across characteristic quartiles",
      compact_row = "Big"
    ),
    lags = lags
  )

  out[["Small_Big_average"]] <- run_ff5_alpha_spread(
    dt = dt,
    y = y_small_big,
    meta = list(
      portfolio_set = portfolio_set,
      construction = construction_name,
      spread_family = spread_family,
      spread_row = "Small-Big",
      spread_definition = "Small average chi spread minus big average chi spread",
      compact_row = "Small-Big"
    ),
    lags = lags
  )

  rbindlist(out, fill = TRUE)
}

chi_32_hl_specs <- list(
  BM_CHI_FIRM = list(
    construction_name = "Firm",
    prefix = chi_32_specs$BM_CHI_FIRM$prefix,
    portfolio_set = "BM_chi",
    spread_family = "32 portfolios: Size-BM-chi"
  ),
  BM_CHI_SIC = list(
    construction_name = "Industry",
    prefix = chi_32_specs$BM_CHI_SIC$prefix,
    portfolio_set = "BM_chi",
    spread_family = "32 portfolios: Size-BM-chi"
  ),
  BM_CHI_SIC2 = list(
    construction_name = "Major group",
    prefix = chi_32_specs$BM_CHI_SIC2$prefix,
    portfolio_set = "BM_chi",
    spread_family = "32 portfolios: Size-BM-chi"
  ),
  OP_CHI_FIRM = list(
    construction_name = "Firm",
    prefix = chi_32_specs$OP_CHI_FIRM$prefix,
    portfolio_set = "OP_chi",
    spread_family = "32 portfolios: Size-OP-chi"
  ),
  OP_CHI_SIC = list(
    construction_name = "Industry",
    prefix = chi_32_specs$OP_CHI_SIC$prefix,
    portfolio_set = "OP_chi",
    spread_family = "32 portfolios: Size-OP-chi"
  ),
  OP_CHI_SIC2 = list(
    construction_name = "Major group",
    prefix = chi_32_specs$OP_CHI_SIC2$prefix,
    portfolio_set = "OP_chi",
    spread_family = "32 portfolios: Size-OP-chi"
  ),
  INV_CHI_FIRM = list(
    construction_name = "Firm",
    prefix = chi_32_specs$INV_CHI_FIRM$prefix,
    portfolio_set = "INV_chi",
    spread_family = "32 portfolios: Size-INV-chi"
  ),
  INV_CHI_SIC = list(
    construction_name = "Industry",
    prefix = chi_32_specs$INV_CHI_SIC$prefix,
    portfolio_set = "INV_chi",
    spread_family = "32 portfolios: Size-INV-chi"
  ),
  INV_CHI_SIC2 = list(
    construction_name = "Major group",
    prefix = chi_32_specs$INV_CHI_SIC2$prefix,
    portfolio_set = "INV_chi",
    spread_family = "32 portfolios: Size-INV-chi"
  )
)

chi_32_hl_ff5_alpha <- rbindlist(
  lapply(chi_32_hl_specs, function(sp) {
    make_32_chi_hl_ff5_alphas(
      dt = factor_test_panel,
      construction_name = sp$construction_name,
      prefix = sp$prefix,
      portfolio_set = sp$portfolio_set,
      spread_family = sp$spread_family,
      lags = NW_LAGS
    )
  }),
  fill = TRUE
)

### 10.42D Combine, format, and save #########################################

chi_hl_ff5_alpha_detailed <- rbindlist(
  list(
    size_chi_hl_ff5_alpha,
    chi_32_hl_ff5_alpha
  ),
  fill = TRUE
)

chi_hl_ff5_alpha_detailed[, alpha_pct := 100 * alpha]

chi_hl_ff5_alpha_detailed[, alpha_display := sprintf(
  "%.2f%s\n(t=%.2f, p=%.3f)",
  alpha_pct,
  alpha_stars,
  alpha_tstat,
  alpha_pval
)]

chi_hl_ff5_alpha_detailed[, alpha_latex := mapply(
  format_alpha_latex_cell,
  alpha = alpha,
  tstat = alpha_tstat,
  stars = alpha_stars
)]

chi_hl_ff5_alpha_compact <- chi_hl_ff5_alpha_detailed[
  !is.na(compact_row),
  .(
    portfolio_set,
    construction,
    spread_family,
    compact_row,
    spread_definition,
    alpha,
    alpha_pct,
    alpha_tstat,
    alpha_pval,
    alpha_stars,
    adj_r2,
    alpha_display,
    alpha_latex
  )
]

portfolio_set_order <- c("Size_chi", "BM_chi", "OP_chi", "INV_chi")
construction_order <- c("Firm", "Industry", "Major group")
compact_row_order <- c("Small", "Big", "Small-Big")

chi_hl_ff5_alpha_compact[, portfolio_set := factor(portfolio_set, levels = portfolio_set_order)]
chi_hl_ff5_alpha_compact[, construction := factor(construction, levels = construction_order)]
chi_hl_ff5_alpha_compact[, compact_row := factor(compact_row, levels = compact_row_order)]

setorder(chi_hl_ff5_alpha_compact, compact_row, portfolio_set, construction)

chi_hl_ff5_alpha_display <- dcast(
  chi_hl_ff5_alpha_compact,
  compact_row ~ portfolio_set + construction,
  value.var = "alpha_latex"
)

setnames(
  chi_hl_ff5_alpha_display,
  old = names(chi_hl_ff5_alpha_display),
  new = gsub(" ", "_", names(chi_hl_ff5_alpha_display))
)

saveRDS(
  chi_hl_ff5_alpha_detailed,
  "section10_chi_hl_ff5_alpha_detailed.rds"
)

saveRDS(
  chi_hl_ff5_alpha_compact,
  "section10_chi_hl_ff5_alpha_compact.rds"
)

saveRDS(
  chi_hl_ff5_alpha_display,
  "section10_chi_hl_ff5_alpha_display.rds"
)

### 10.42E Console output #####################################################

cat("\n============================================================\n")
cat("FF5 ALPHAS OF HIGH-MINUS-LOW CHI SPREAD PORTFOLIOS\n")
cat("============================================================\n\n")

cat("Compact numeric table:\n")
print(round_dt(chi_hl_ff5_alpha_compact, 4))

cat("\nLaTeX-ready display table:\n")
print(chi_hl_ff5_alpha_display)





### 10.43 LaTeX Tables for Chi Alpha Results ##################################
### 10.43A Shared formatter ###################################################

format_alpha_latex_cell <- function(alpha, tstat, stars, digits_est = 2, digits_t = 2, scale = 100) {
  if (is.na(alpha) | is.na(tstat)) {
    return("")
  }

  star_latex <- ifelse(
    is.na(stars) | stars == "",
    "",
    paste0("^{", stars, "}")
  )

  sprintf(
    paste0("\\makecell{$%.", digits_est, "f%s$ \\\\ $(%.", digits_t, "f)$}"),
    scale * alpha,
    star_latex,
    tstat
  )
}

### 10.43B Table: FF5 alphas across 25 Size-chi portfolios ####################

size_chi_ff5_alpha_profile <- copy(
  model_alpha_chi_25_vw[
    model == "FF5" &
      portfolio_set %in% c("CHI_FIRM_25", "CHI_SIC_25", "CHI_SIC2_25")
  ]
)

size_chi_ff5_alpha_profile[, construction := fcase(
  portfolio_set == "CHI_FIRM_25", "Firm",
  portfolio_set == "CHI_SIC_25", "Industry",
  portfolio_set == "CHI_SIC2_25", "Major group",
  default = NA_character_
)]

size_chi_ff5_alpha_profile[, size_bucket := factor(size_bucket, levels = paste0("S", 1:5))]
size_chi_ff5_alpha_profile[, chi_bucket := factor(chi_bucket, levels = paste0("C", 1:5))]
size_chi_ff5_alpha_profile[, construction := factor(construction, levels = c("Firm", "Industry", "Major group"))]

size_chi_ff5_alpha_profile[, alpha_latex := mapply(
  format_alpha_latex_cell,
  alpha = alpha,
  tstat = alpha_tstat,
  stars = alpha_stars
)]

write_size_chi_ff5_alpha_table <- function(
    dt,
    file = "table_size_chi_ff5_alpha_profile.tex") {
  con_order <- c("Firm", "Industry", "Major group")
  con_titles <- c(
    Firm = "Firm-level $\\chi$",
    Industry = "Industry $\\chi$",
    `Major group` = "Major-group $\\chi$"
  )

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Table: FF5 alphas across Size-chi portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\caption{FF5 alphas across Size--$\\chi$ portfolios}\n")
  cat("\\label{tab:size_chi_ff5_alpha_profile}\n\n")

  cat("\\begin{minipage}{0.90\\linewidth}\n")
  cat("\\footnotesize\n")
  cat("US equities, July 1982--December 2024. This table reports FF5 alphas for the 25 Size--$\\chi$ portfolios. Rows correspond to size quintiles, where S1 denotes the smallest firms and S5 denotes the largest firms. Columns correspond to $\\chi$ quintiles, where C1 denotes the lowest-$\\chi$ portfolios and C5 denotes the highest-$\\chi$ portfolios. Panel A uses firm-level $\\chi$, Panel B uses SIC 4-digit industry-level $\\chi$, and Panel C uses SIC 2-digit major-group $\\chi$. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{minipage}\n\n")

  cat("\\vspace{0.7em}\n\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{8pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.20}\n\n")

  cat("\\begin{tabular}{lccccc}\n")
  cat("\\toprule\n")
  cat("Size & C1 & C2 & C3 & C4 & C5 \\\\\n")
  cat("\\midrule\n")

  for (kk in seq_along(con_order)) {
    con <- con_order[kk]
    panel_letter <- LETTERS[kk]

    cat("\\multicolumn{6}{l}{\\textit{Panel ", panel_letter, ": ", con_titles[con], "}} \\\\\n", sep = "")
    cat("\\midrule\n")

    for (ss in paste0("S", 1:5)) {
      cells <- sapply(paste0("C", 1:5), function(cc) {
        val <- dt[
          construction == con &
            size_bucket == ss &
            chi_bucket == cc,
          alpha_latex
        ]
        if (length(val) == 0) "" else val[1]
      })

      cat(ss, " & ", paste(cells, collapse = " & "), " \\\\\n", sep = "")
    }

    if (kk < length(con_order)) {
      cat("\\addlinespace\n")
      cat("\\midrule\n")
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")

  sink()
}

write_size_chi_ff5_alpha_table(size_chi_ff5_alpha_profile)

saveRDS(
  size_chi_ff5_alpha_profile,
  "section10_size_chi_ff5_alpha_profile.rds"
)

### 10.43C Table: FF5 alphas of high-minus-low chi spreads ####################

write_chi_hl_ff5_alpha_spread_table <- function(
    dt,
    file = "table_chi_hl_ff5_alpha_spreads.tex") {
  x <- copy(dt)
  x[, portfolio_set := as.character(portfolio_set)]
  x[, construction := as.character(construction)]
  x[, compact_row := as.character(compact_row)]

  get_cell <- function(row_name, port_name, con_name) {
    val <- x[
      compact_row == row_name &
        portfolio_set == port_name &
        construction == con_name,
      alpha_latex
    ]
    if (length(val) == 0) "" else val[1]
  }

  make_row <- function(row_name, row_label) {
    cells <- c(
      get_cell(row_name, "Size_chi", "Firm"),
      get_cell(row_name, "Size_chi", "Industry"),
      get_cell(row_name, "Size_chi", "Major group"),
      get_cell(row_name, "BM_chi", "Firm"),
      get_cell(row_name, "BM_chi", "Industry"),
      get_cell(row_name, "BM_chi", "Major group"),
      get_cell(row_name, "OP_chi", "Firm"),
      get_cell(row_name, "OP_chi", "Industry"),
      get_cell(row_name, "OP_chi", "Major group"),
      get_cell(row_name, "INV_chi", "Firm"),
      get_cell(row_name, "INV_chi", "Industry"),
      get_cell(row_name, "INV_chi", "Major group")
    )

    paste0(row_label, " & ", paste(cells, collapse = " & "), " \\\\")
  }

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Table: FF5 alphas of high-minus-low chi spread portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{landscape}\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\caption{FF5 alphas of high-minus-low $\\chi$ spread portfolios}\n")
  cat("\\label{tab:chi_hl_ff5_alpha_spreads}\n\n")

  cat("\\begin{minipage}{0.97\\linewidth}\n")
  cat("\\footnotesize\n")
  cat("US equities, July 1982--December 2024. This table reports FF5 alphas from high-minus-low $\\chi$ spread portfolios. For the 25 Size--$\\chi$ portfolios, Small denotes the highest-minus-lowest $\\chi$ spread within the smallest size quintile, while Big denotes the corresponding spread within the largest size quintile. For the 32-portfolio families, Small and Big denote average highest-minus-lowest $\\chi$ spreads across the four non-$\\chi$ characteristic quartiles. Small--Big is the difference between the small-firm and big-firm $\\chi$ spreads. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{minipage}\n\n")

  cat("\\vspace{0.7em}\n\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{3pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.20}\n\n")

  cat("\\begin{tabular*}{0.97\\linewidth}{@{\\extracolsep{\\fill}}lcccccccccccc@{}}\n")
  cat("\\toprule\n")
  cat("& \\multicolumn{3}{c}{25 portfolios: Size--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--B/M--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--OP--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--INV--$\\chi$} \\\\\n")
  cat("\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group \\\\\n")
  cat("\\midrule\n")

  cat(make_row("Small", "Small"), "\n")
  cat(make_row("Big", "Big"), "\n")
  cat(make_row("Small-Big", "Small--Big"), "\n")

  cat("\\bottomrule\n")
  cat("\\end{tabular*}\n")
  cat("\\end{table}\n")
  cat("\\end{landscape}\n")

  sink()
}

write_chi_hl_ff5_alpha_spread_table(chi_hl_ff5_alpha_compact)

saveRDS(
  chi_hl_ff5_alpha_compact,
  "section10_chi_hl_ff5_alpha_spread_table_input.rds"
)

### 10.43D Console confirmation ###############################################

cat("\n============================================================\n")
cat("LATEX TABLES CREATED\n")
cat("============================================================\n\n")
cat("Created files:\n")
cat("- table_size_chi_ff5_alpha_profile.tex\n")
cat("- table_chi_hl_ff5_alpha_spreads.tex\n\n")


size_chi_c1_c5_check <- size_chi_ff5_alpha_profile[
  ,
  .(
    C1_alpha = alpha[chi_bucket == "C1"],
    C5_alpha = alpha[chi_bucket == "C5"],
    C5_minus_C1 = alpha[chi_bucket == "C5"] - alpha[chi_bucket == "C1"],
    C1_t = alpha_tstat[chi_bucket == "C1"],
    C5_t = alpha_tstat[chi_bucket == "C5"]
  ),
  by = .(construction, size_bucket)
]

print(round_dt(size_chi_c1_c5_check, 4))


### 10.44 Display Generated LaTeX Tables #####################################
### Check that the table files exist
table_files <- c(
  "table_size_chi_ff5_alpha_profile.tex",
  "table_chi_hl_ff5_alpha_spreads.tex"
)

print(file.exists(table_files))
print(table_files)


### Display the Size-chi FF5 alpha profile table
cat("\n\n")
cat("TABLE: FF5 alphas across Size-chi portfolios\n")
cat("\n\n")

cat(
  readLines("table_size_chi_ff5_alpha_profile.tex"),
  sep = "\n"
)


### Display the high-minus-low chi spread alpha table
cat("\n\n")
cat("TABLE: FF5 alphas of high-minus-low chi spread portfolios\n")
cat("\n\n")

cat(
  readLines("table_chi_hl_ff5_alpha_spreads.tex"),
  sep = "\n"
)


### Open the generated .tex files in RStudio viewer/editor
file.show("table_size_chi_ff5_alpha_profile.tex")
file.show("table_chi_hl_ff5_alpha_spreads.tex")




### 10.44 Robustness Tables for Chi Alpha Results #############################
### 10.44A Shared LaTeX formatter #############################################

format_alpha_latex_cell <- function(alpha, tstat, stars, digits_est = 2, digits_t = 2, scale = 100) {
  if (is.na(alpha) | is.na(tstat)) {
    return("")
  }

  star_latex <- ifelse(
    is.na(stars) | stars == "",
    "",
    paste0("^{", stars, "}")
  )

  sprintf(
    paste0("\\makecell{$%.", digits_est, "f%s$ \\\\ $(%.", digits_t, "f)$}"),
    scale * alpha,
    star_latex,
    tstat
  )
}

### 10.44B Robustness Table A: Size-chi FF5 alpha profile #####################

write_robustness_size_chi_alpha_profile <- function(
    dt,
    file = "table_robustness_size_chi_ff5_alpha_profile.tex") {
  con_order <- c("Firm", "Industry", "Major group")
  con_titles <- c(
    Firm = "Firm-level $\\chi$",
    Industry = "Industry $\\chi$",
    `Major group` = "Major-group $\\chi$"
  )

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Robustness Table: FF5 alphas across Size-chi portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\caption{Robustness: FF5 alphas across Size--$\\chi$ portfolios}\n")
  cat("\\label{tab:robustness_size_chi_ff5_alpha_profile}\n\n")

  cat("\\begin{minipage}{0.90\\linewidth}\n")
  cat("\\footnotesize\n")
  cat("US equities, July 1982--December 2024. This table reports FF5 alphas for the 25 Size--$\\chi$ portfolios across alternative estimates of $\\chi$. Rows correspond to size quintiles, where S1 denotes the smallest firms and S5 denotes the largest firms. Columns correspond to $\\chi$ quintiles, where C1 denotes the lowest-$\\chi$ portfolios and C5 denotes the highest-$\\chi$ portfolios. Panel A uses firm-level $\\chi$, Panel B uses SIC 4-digit industry-level $\\chi$, and Panel C uses SIC 2-digit major-group $\\chi$. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{minipage}\n\n")

  cat("\\vspace{0.7em}\n\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{8pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.20}\n\n")
  cat("\\begin{tabular}{lccccc}\n")
  cat("\\toprule\n")
  cat("Size & C1 & C2 & C3 & C4 & C5 \\\\\n")
  cat("\\midrule\n")

  for (kk in seq_along(con_order)) {
    con <- con_order[kk]
    panel_letter <- LETTERS[kk]

    cat("\\multicolumn{6}{l}{\\textit{Panel ", panel_letter, ": ", con_titles[con], "}} \\\\\n", sep = "")
    cat("\\midrule\n")

    for (ss in paste0("S", 1:5)) {
      cells <- sapply(paste0("C", 1:5), function(cc) {
        val <- dt[
          construction == con &
            size_bucket == ss &
            chi_bucket == cc,
          alpha_latex
        ]
        if (length(val) == 0) "" else val[1]
      })

      cat(ss, " & ", paste(cells, collapse = " & "), " \\\\\n", sep = "")
    }

    if (kk < length(con_order)) {
      cat("\\addlinespace\n")
      cat("\\midrule\n")
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")

  sink()
}

write_robustness_size_chi_alpha_profile(size_chi_ff5_alpha_profile)

### 10.44C Robustness Table B: H-L chi FF5 alpha spreads ######################

write_robustness_chi_hl_alpha_spreads <- function(
    dt,
    file = "table_robustness_chi_hl_ff5_alpha_spreads.tex") {
  x <- copy(dt)
  x[, portfolio_set := as.character(portfolio_set)]
  x[, construction := as.character(construction)]
  x[, compact_row := as.character(compact_row)]

  get_cell <- function(row_name, port_name, con_name) {
    val <- x[
      compact_row == row_name &
        portfolio_set == port_name &
        construction == con_name,
      alpha_latex
    ]
    if (length(val) == 0) "" else val[1]
  }

  make_row <- function(row_name, row_label) {
    cells <- c(
      get_cell(row_name, "Size_chi", "Firm"),
      get_cell(row_name, "Size_chi", "Industry"),
      get_cell(row_name, "Size_chi", "Major group"),
      get_cell(row_name, "BM_chi", "Firm"),
      get_cell(row_name, "BM_chi", "Industry"),
      get_cell(row_name, "BM_chi", "Major group"),
      get_cell(row_name, "OP_chi", "Firm"),
      get_cell(row_name, "OP_chi", "Industry"),
      get_cell(row_name, "OP_chi", "Major group"),
      get_cell(row_name, "INV_chi", "Firm"),
      get_cell(row_name, "INV_chi", "Industry"),
      get_cell(row_name, "INV_chi", "Major group")
    )

    paste0(row_label, " & ", paste(cells, collapse = " & "), " \\\\")
  }

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Robustness Table: FF5 alphas of H-L chi spread portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{landscape}\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\caption{Robustness: FF5 alphas of high-minus-low $\\chi$ spread portfolios}\n")
  cat("\\label{tab:robustness_chi_hl_ff5_alpha_spreads}\n\n")

  cat("\\begin{minipage}{0.97\\linewidth}\n")
  cat("\\footnotesize\n")
  cat("US equities, July 1982--December 2024. This table reports FF5 alphas from high-minus-low $\\chi$ spread portfolios across alternative estimates of $\\chi$. For the 25 Size--$\\chi$ portfolios, Small denotes the highest-minus-lowest $\\chi$ spread within the smallest size quintile, while Big denotes the corresponding spread within the largest size quintile. For the 32-portfolio families, Small and Big denote average highest-minus-lowest $\\chi$ spreads across the four non-$\\chi$ characteristic quartiles. Small--Big is the difference between the small-firm and big-firm $\\chi$ spreads. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{minipage}\n\n")

  cat("\\vspace{0.7em}\n\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{3pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.20}\n\n")
  cat("\\begin{tabular*}{0.97\\linewidth}{@{\\extracolsep{\\fill}}lcccccccccccc@{}}\n")
  cat("\\toprule\n")
  cat("& \\multicolumn{3}{c}{25 portfolios: Size--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--B/M--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--OP--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{32 portfolios: Size--INV--$\\chi$} \\\\\n")
  cat("\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group \\\\\n")
  cat("\\midrule\n")

  cat(make_row("Small", "Small"), "\n")
  cat(make_row("Big", "Big"), "\n")
  cat(make_row("Small-Big", "Small--Big"), "\n")

  cat("\\bottomrule\n")
  cat("\\end{tabular*}\n")
  cat("\\end{table}\n")
  cat("\\end{landscape}\n")

  sink()
}

write_robustness_chi_hl_alpha_spreads(chi_hl_ff5_alpha_compact)

### 10.44D Optional summary table for writing interpretation ##################

robustness_chi_hl_summary <- chi_hl_ff5_alpha_compact[
  ,
  .(
    portfolio_set,
    construction,
    compact_row,
    alpha_pct,
    alpha_tstat,
    alpha_pval,
    alpha_stars,
    adj_r2
  )
][order(compact_row, portfolio_set, construction)]

saveRDS(
  robustness_chi_hl_summary,
  "section10_robustness_chi_hl_summary.rds"
)

### 10.44E Console confirmation ###############################################

cat("\n============================================================\n")
cat("ROBUSTNESS TABLES CREATED\n")
cat("============================================================\n\n")

cat("Created files:\n")
cat("- table_robustness_size_chi_ff5_alpha_profile.tex\n")
cat("- table_robustness_chi_hl_ff5_alpha_spreads.tex\n\n")

cat("Main robustness summary:\n")
print(round_dt(robustness_chi_hl_summary, 4))


### 10.45 Pretty LaTeX Robustness Tables ######################################

### Required LaTeX packages:
### \usepackage{booktabs}
### \usepackage{makecell}
### \usepackage{threeparttable}
### \usepackage{pdflscape}
### \usepackage{adjustbox}

### 10.45A Pretty Size-chi FF5 alpha profile table ############################

write_pretty_size_chi_alpha_profile <- function(
    dt,
    file = "table_pretty_size_chi_ff5_alpha_profile.tex") {
  con_order <- c("Firm", "Industry", "Major group")
  con_titles <- c(
    Firm = "Firm-level $\\chi$",
    Industry = "Industry-level $\\chi$",
    `Major group` = "Major-group $\\chi$"
  )

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Pretty Table: FF5 alphas across Size-chi portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\begin{threeparttable}\n")
  cat("\\caption{FF5 alphas across Size--$\\chi$ portfolios}\n")
  cat("\\label{tab:size_chi_ff5_alpha_profile}\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{7pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.18}\n\n")

  cat("\\begin{tabular}{lccccc}\n")
  cat("\\toprule\n")
  cat("& \\multicolumn{5}{c}{$\\chi$ quintile} \\\\\n")
  cat("\\cmidrule(lr){2-6}\n")
  cat("Size & C1 & C2 & C3 & C4 & C5 \\\\\n")
  cat("\\midrule\n")

  for (kk in seq_along(con_order)) {
    con <- con_order[kk]
    panel_letter <- LETTERS[kk]

    cat("\\multicolumn{6}{l}{\\textit{Panel ", panel_letter, ". ", con_titles[con], "}} \\\\\n", sep = "")
    cat("\\addlinespace[0.15em]\n")

    for (ss in paste0("S", 1:5)) {
      cells <- sapply(paste0("C", 1:5), function(cc) {
        val <- dt[
          construction == con &
            size_bucket == ss &
            chi_bucket == cc,
          alpha_latex
        ]
        if (length(val) == 0) "" else val[1]
      })

      cat(ss, " & ", paste(cells, collapse = " & "), " \\\\\n", sep = "")
    }

    if (kk < length(con_order)) {
      cat("\\addlinespace[0.45em]\n")
    }
  }

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n\n")

  cat("\\begin{tablenotes}[flushleft]\n")
  cat("\\footnotesize\n")
  cat("\\item \\textit{Notes:} This table reports FF5 alphas for the 25 Size--$\\chi$ portfolios. Rows correspond to size quintiles, where S1 denotes the smallest firms and S5 denotes the largest firms. Columns correspond to $\\chi$ quintiles, where C1 denotes the lowest-$\\chi$ portfolios and C5 denotes the highest-$\\chi$ portfolios. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat("\\end{table}\n")

  sink()
}

write_pretty_size_chi_alpha_profile(size_chi_ff5_alpha_profile)


### 10.45B Pretty high-minus-low chi spread alpha table #######################

write_pretty_chi_hl_alpha_spreads <- function(
    dt,
    file = "table_pretty_chi_hl_ff5_alpha_spreads.tex") {
  x <- copy(dt)
  x[, portfolio_set := as.character(portfolio_set)]
  x[, construction := as.character(construction)]
  x[, compact_row := as.character(compact_row)]

  get_cell <- function(row_name, port_name, con_name) {
    val <- x[
      compact_row == row_name &
        portfolio_set == port_name &
        construction == con_name,
      alpha_latex
    ]
    if (length(val) == 0) "" else val[1]
  }

  make_row <- function(row_name, row_label) {
    cells <- c(
      get_cell(row_name, "Size_chi", "Firm"),
      get_cell(row_name, "Size_chi", "Industry"),
      get_cell(row_name, "Size_chi", "Major group"),
      get_cell(row_name, "BM_chi", "Firm"),
      get_cell(row_name, "BM_chi", "Industry"),
      get_cell(row_name, "BM_chi", "Major group"),
      get_cell(row_name, "OP_chi", "Firm"),
      get_cell(row_name, "OP_chi", "Industry"),
      get_cell(row_name, "OP_chi", "Major group"),
      get_cell(row_name, "INV_chi", "Firm"),
      get_cell(row_name, "INV_chi", "Industry"),
      get_cell(row_name, "INV_chi", "Major group")
    )

    paste0("\\textbf{", row_label, "} & ", paste(cells, collapse = " & "), " \\\\")
  }

  sink(file)

  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("% Pretty Table: FF5 alphas of high-minus-low chi spread portfolios\n")
  cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  cat("\\begin{landscape}\n")
  cat("\\begin{table}[p]\n")
  cat("\\centering\n")
  cat("\\begin{threeparttable}\n")
  cat("\\caption{FF5 alphas of high-minus-low $\\chi$ spread portfolios}\n")
  cat("\\label{tab:chi_hl_ff5_alpha_spreads}\n")
  cat("\\scriptsize\n")
  cat("\\setlength{\\tabcolsep}{3pt}\n")
  cat("\\renewcommand{\\arraystretch}{1.20}\n\n")

  cat("\\begin{adjustbox}{max width=\\linewidth}\n")
  cat("\\begin{tabular}{lcccccccccccc}\n")
  cat("\\toprule\n")
  cat("& \\multicolumn{3}{c}{Size--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{Size--B/M--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{Size--OP--$\\chi$} \n")
  cat("& \\multicolumn{3}{c}{Size--INV--$\\chi$} \\\\\n")
  cat("\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10} \\cmidrule(lr){11-13}\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group\n")
  cat("& Firm & Industry & Major group \\\\\n")
  cat("\\midrule\n")

  cat(make_row("Small", "Small"), "\n")
  cat("\\addlinespace[0.25em]\n")
  cat(make_row("Big", "Big"), "\n")
  cat("\\addlinespace[0.25em]\n")
  cat(make_row("Small-Big", "Small--Big"), "\n")

  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{adjustbox}\n\n")

  cat("\\begin{tablenotes}[flushleft]\n")
  cat("\\footnotesize\n")
  cat("\\item \\textit{Notes:} This table reports FF5 alphas from high-minus-low $\\chi$ spread portfolios. For the 25 Size--$\\chi$ portfolios, Small denotes the highest-minus-lowest $\\chi$ spread within the smallest size quintile, while Big denotes the corresponding spread within the largest size quintile. For the 32-portfolio families, Small and Big denote average highest-minus-lowest $\\chi$ spreads across the four non-$\\chi$ characteristic quartiles. Small--Big is the difference between the small-firm and big-firm $\\chi$ spreads. Alphas are reported in monthly percent. Newey--West $t$-statistics with 12 lags are reported in parentheses. Statistical significance is denoted by $^{*}$, $^{**}$, and $^{***}$ for the 10\\%, 5\\%, and 1\\% levels, respectively.\n")
  cat("\\end{tablenotes}\n")
  cat("\\end{threeparttable}\n")
  cat("\\end{table}\n")
  cat("\\end{landscape}\n")

  sink()
}

write_pretty_chi_hl_alpha_spreads(chi_hl_ff5_alpha_compact)



### 10.46 H-L chi alpha spreads across alternative models #####################

### 10.46A Model helper #######################################################

get_hl_model_rhs <- function(model_name, variant_name) {
  if (model_name == "CAPM") {
    return(c("MKT"))
  }

  if (model_name == "FF3") {
    return(c("MKT", "SMB", "HML"))
  }

  if (model_name == "FF5") {
    return(c("MKT", "SMB", "HML", "RMW", "CMA"))
  }

  if (model_name == "FF5_SMQ") {
    return(variant_factor_map[[variant_name]]$ff5_smq)
  }

  stop("Unknown model: ", model_name)
}

### 10.46B Run H-L regressions ################################################

run_size_chi_hl_models <- function(
    dt,
    construction_name,
    prefix,
    variant_name,
    model_names = c("CAPM", "FF3", "FF5", "FF5_SMQ"),
    lags = NW_LAGS) {
  out <- list()

  for (model_name in model_names) {
    rhs <- get_hl_model_rhs(model_name, variant_name)

    for (ss in 1:5) {
      low_col <- paste0(prefix, "S", ss, "C1")
      high_col <- paste0(prefix, "S", ss, "C5")

      if (!low_col %in% names(dt)) {
        stop("Missing low-chi column: ", low_col)
      }

      if (!high_col %in% names(dt)) {
        stop("Missing high-chi column: ", high_col)
      }

      y <- dt[[high_col]] - dt[[low_col]]

      reg <- run_ts_reg_nw(
        y = y,
        X = dt[, ..rhs],
        lags = lags
      )

      if (is.null(reg)) {
        next
      }

      a <- reg[term == "(Intercept)"][1]

      out[[paste(construction_name, model_name, ss, sep = "_")]] <-
        data.table(
          construction = construction_name,
          variant = variant_name,
          model = model_name,
          size_bucket = paste0("S", ss),
          spread = "C5-C1",
          alpha = a$estimate,
          alpha_pct = 100 * a$estimate,
          tstat = a$tstat,
          pval = a$pval,
          stars = a$stars,
          adj_r2 = a$adj_r2
        )
    }
  }

  rbindlist(out, fill = TRUE)
}

### 10.46C Apply to firm, industry, and major-group chi ########################

size_chi_hl_model_specs <- list(
  Firm = list(
    prefix = "CHI_FIRM_",
    variant_name = "FIRM_VW"
  ),
  Industry = list(
    prefix = "CHI_SIC_",
    variant_name = "SIC_VW"
  ),
  `Major group` = list(
    prefix = "CHI_SIC2_",
    variant_name = "SIC2_VW"
  )
)

size_chi_hl_alternative_models <- rbindlist(
  lapply(names(size_chi_hl_model_specs), function(nm) {
    sp <- size_chi_hl_model_specs[[nm]]

    run_size_chi_hl_models(
      dt = factor_test_panel,
      construction_name = nm,
      prefix = sp$prefix,
      variant_name = sp$variant_name,
      model_names = c("CAPM", "FF3", "FF5", "FF5_SMQ"),
      lags = NW_LAGS
    )
  }),
  fill = TRUE
)

### 10.46D Display and save ###################################################

size_chi_hl_alternative_models[, display := sprintf(
  "%.2f%s\n(%.2f)",
  alpha_pct,
  stars,
  tstat
)]

size_chi_hl_alternative_models_display <- dcast(
  size_chi_hl_alternative_models,
  construction + size_bucket ~ model,
  value.var = "display"
)

saveRDS(
  size_chi_hl_alternative_models,
  "section10_size_chi_hl_alternative_models.rds"
)

saveRDS(
  size_chi_hl_alternative_models_display,
  "section10_size_chi_hl_alternative_models_display.rds"
)

cat("\n============================================================\n")
cat("SIZE-CHI H-L ALPHAS ACROSS ALTERNATIVE MODELS CREATED\n")
cat("============================================================\n\n")

cat("Detailed numeric table:\n")
print(round_dt(size_chi_hl_alternative_models, 4))

cat("\nDisplay table:\n")
print(size_chi_hl_alternative_models_display)




### 10.47 H-L chi FF5 alphas for 32 portfolios ################################
### Matches section 8 specification exactly ####################################

### 10.47A Helper ##############################################################

pretty_char_label <- function(char_name, quartile) {
  labels <- c("Low", "2", "3", "High")
  labels[quartile]
}

run_32_hl_alpha_ff5 <- function(sort_family,
                                construction,
                                char_name,
                                prefix,
                                lags = NW_LAGS) {
  rhs <- c("MKT", "SMB", "HML", "RMW", "CMA")

  out <- list()

  for (size_code in c("S", "B")) {
    for (q in 1:4) {
      low_col <- paste0(prefix, size_code, q, "1")
      high_col <- paste0(prefix, size_code, q, "4")

      if (!low_col %in% names(factor_test_panel)) {
        stop("Missing low-chi column: ", low_col)
      }

      if (!high_col %in% names(factor_test_panel)) {
        stop("Missing high-chi column: ", high_col)
      }

      y_hl <- factor_test_panel[[high_col]] - factor_test_panel[[low_col]]

      reg <- run_ts_reg_nw(
        y    = y_hl,
        X    = factor_test_panel[, ..rhs],
        lags = lags
      )

      if (is.null(reg)) next

      a <- reg[term == "(Intercept)"][1]

      out[[paste(sort_family, construction, size_code, q, sep = "_")]] <-
        data.table(
          sort_family       = sort_family,
          construction      = construction,
          model             = "FF5",
          char_name         = char_name,
          size_code         = size_code,
          size_bucket       = fifelse(size_code == "S", "Small", "Big"),
          char_bucket       = paste0("Q", q),
          char_bucket_label = pretty_char_label(char_name, q),
          spread            = "Q4-Q1 chi",
          low_col           = low_col,
          high_col          = high_col,
          alpha             = a$estimate,
          alpha_pct         = 100 * a$estimate,
          tstat             = a$tstat,
          pval              = a$pval,
          stars             = a$stars,
          adj_r2            = a$adj_r2
        )
    }
  }

  rbindlist(out, fill = TRUE)
}

### 10.47B Prefix specs ########################################################

sort_family_order <- c("Size-OP-chi", "Size-BM-chi", "Size-INV-chi")
sort_order_dt <- data.table(
  sort_family = sort_family_order,
  sort_order = seq_along(sort_family_order)
)
size_order_dt <- data.table(
  size_bucket = c("Small", "Big"),
  size_order = 1:2
)
char_order_dt <- data.table(
  char_bucket = paste0("Q", 1:4),
  char_order = 1:4
)

specs_1047 <- CJ(
  sort_family = sort_family_order,
  construction = construction_order,
  unique = TRUE
)

specs_1047[, char_name := fcase(
  sort_family == "Size-OP-chi",  "OP",
  sort_family == "Size-BM-chi",  "B/M",
  sort_family == "Size-INV-chi", "Inv"
)]

specs_1047[, prefix := fcase(
  sort_family == "Size-OP-chi" & construction == "Firm", "OP_CHI_FIRM_",
  sort_family == "Size-OP-chi" & construction == "Industry", "OP_CHI_SIC_",
  sort_family == "Size-OP-chi" & construction == "Major group", "OP_CHI_SIC2_",
  sort_family == "Size-BM-chi" & construction == "Firm", "BM_CHI_FIRM_",
  sort_family == "Size-BM-chi" & construction == "Industry", "BM_CHI_SIC_2X4X4_",
  sort_family == "Size-BM-chi" & construction == "Major group", "BM_CHI_SIC2_",
  sort_family == "Size-INV-chi" & construction == "Firm", "INV_CHI_FIRM_",
  sort_family == "Size-INV-chi" & construction == "Industry", "INV_CHI_SIC_2X4X4_",
  sort_family == "Size-INV-chi" & construction == "Major group", "INV_CHI_SIC2_"
)]

cat("\n============================================================\n")
cat("10.47 PREFIXES USED\n")
cat("============================================================\n")
print(specs_1047)

### 10.47C Run regressions #####################################################

chi_32_hl_ff5_alphas <- rbindlist(
  lapply(seq_len(nrow(specs_1047)), function(k) {
    run_32_hl_alpha_ff5(
      sort_family  = specs_1047$sort_family[k],
      construction = specs_1047$construction[k],
      char_name    = specs_1047$char_name[k],
      prefix       = specs_1047$prefix[k],
      lags         = NW_LAGS
    )
  }),
  fill = TRUE
)

### 10.47D Display objects #####################################################

chi_32_hl_ff5_alphas[, display := sprintf(
  "%.2f%s\n(%.2f)",
  alpha_pct,
  stars,
  tstat
)]

latex_stars <- function(x) {
  fifelse(x == "", "", paste0("$^{", x, "}$"))
}

chi_32_hl_ff5_alphas[, display_latex := sprintf(
  "\\makecell{%.2f%s \\\\ (%.2f)}",
  alpha_pct,
  latex_stars(stars),
  tstat
)]

chi_32_hl_ff5_display <- dcast(
  chi_32_hl_ff5_alphas,
  sort_family + char_name + size_bucket + char_bucket + char_bucket_label ~ construction,
  value.var = "display"
)

chi_32_hl_ff5_latex_display <- dcast(
  chi_32_hl_ff5_alphas,
  sort_family + char_name + size_bucket + char_bucket + char_bucket_label ~ construction,
  value.var = "display_latex"
)

### Sort both display tables ###################################################

for (disp_dt in list(chi_32_hl_ff5_display, chi_32_hl_ff5_latex_display)) {
  disp_dt <- merge(disp_dt, sort_order_dt, by = "sort_family")
  disp_dt <- merge(disp_dt, size_order_dt, by = "size_bucket")
  disp_dt <- merge(disp_dt, char_order_dt, by = "char_bucket")
  setorder(disp_dt, sort_order, char_order, size_order)
  disp_dt[, c("sort_order", "size_order", "char_order") := NULL]
}

setcolorder(
  chi_32_hl_ff5_display,
  c(
    "sort_family", "char_name", "size_bucket", "char_bucket",
    "char_bucket_label", construction_order
  )
)

setcolorder(
  chi_32_hl_ff5_latex_display,
  c(
    "sort_family", "char_name", "size_bucket", "char_bucket",
    "char_bucket_label", construction_order
  )
)

### 10.47E Save ################################################################

saveRDS(
  chi_32_hl_ff5_alphas,
  "section10_47_chi_32_hl_ff5_alphas.rds"
)

saveRDS(
  chi_32_hl_ff5_display,
  "section10_47_chi_32_hl_ff5_display.rds"
)

saveRDS(
  chi_32_hl_ff5_latex_display,
  "section10_47_chi_32_hl_ff5_latex_display.rds"
)

### 10.47F Console output ######################################################

cat("\n============================================================\n")
cat("10.47 H-L CHI FF5 ALPHAS FOR 32 PORTFOLIOS\n")
cat("============================================================\n\n")

cat("Detailed numeric table:\n")
print(round_dt(chi_32_hl_ff5_alphas, 4))

cat("\nDisplay table:\n")
print(chi_32_hl_ff5_display)
