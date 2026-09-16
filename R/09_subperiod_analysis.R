## 9. Subperiod Evidence on the Investment-q Relationship #####################
## Baseline specification: SIC value-weighted SMQ only #########################

cat("\014")
rm(list = ls())
graphics.off()

### 9.1 Setup #################################################################
library(data.table)
library(lubridate)
library(zoo)
library(sandwich)
library(lmtest)

source(file.path("R", "00_config.R"))

START_DATE <- as.Date("1963-07-01")
END_DATE <- as.Date("2024-12-01")
NW_LAGS <- 12
ROLL_WIN <- 60

### 9.2 Load Data #############################################################
data <- readRDS("crsp_compustat_with_adj_cost.rds")
ff5 <- readRDS("ff5_factors_internal.rds")
SMQ <- readRDS("SMQ_factors.rds")

setDT(data)
setDT(ff5)
setDT(SMQ)

data[, mdate := as.Date(mdate)]
ff5[, mdate := as.Date(mdate)]
SMQ[, mdate := as.Date(mdate)]

### 9.3 Restrict Sample #######################################################
data <- data[
  mdate >= START_DATE &
    mdate <= END_DATE
]

ff5 <- ff5[
  mdate >= START_DATE &
    mdate <= END_DATE
]

SMQ <- SMQ[
  mdate >= START_DATE &
    mdate <= END_DATE
]

### 9.4 Calendar Fields #######################################################
data[, year := year(mdate)]
data[, month := month(mdate)]
data[, ffyear := fifelse(month >= 7, year, year - 1)]

### 9.5 Subperiod Definitions #################################################
subperiods <- list(
  Early = list(
    start = as.Date("1982-07-01"),
    end   = as.Date("2002-12-01")
  ),
  Late = list(
    start = as.Date("2003-01-01"),
    end   = as.Date("2024-12-01")
  )
)

### 9.6 Helpers ###############################################################
ff_30_70_breaks <- function(x) {
  qs <- quantile(x, probs = c(0.3, 0.7), type = 1, na.rm = TRUE)
  c(low = as.numeric(qs[1]), high = as.numeric(qs[2]))
}

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

round_dt <- function(dt, digits = 4) {
  out <- copy(dt)
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  out
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

run_ts_reg_nw <- function(formula, dt, lags = NW_LAGS) {
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

  out <- data.table(
    term     = rownames(ct),
    estimate = as.numeric(ct[, 1]),
    tstat    = as.numeric(ct[, 3]),
    pval     = as.numeric(ct[, 4]),
    stars    = star_p(as.numeric(ct[, 4]))
  )
  out[, adj_r2 := summary(fit)$adj.r.squared]
  out
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

### 9.7 FF-Style Annual Signals ###############################################
op_ff <- unique(data[, .(gvkey, fyear, revt, cogs, xsga, xint, be_clean)])
op_ff[, op_ff := fifelse(
  !is.na(be_clean) & be_clean > 0 & !is.na(revt) & !is.na(cogs),
  (revt - cogs - fifelse(is.na(xsga), 0, xsga) - fifelse(is.na(xint), 0, xint)) / be_clean,
  NA_real_
)]
op_ff[, fyear := fyear + 1L]

inv_ff <- unique(data[, .(gvkey, fyear, at)])
setorder(inv_ff, gvkey, fyear)
inv_ff[, at_lag := shift(at), by = gvkey]
inv_ff[, inv_ff := fifelse(
  !is.na(at_lag) & at_lag > 0,
  (at - at_lag) / at_lag,
  NA_real_
)]
inv_ff[, fyear := fyear + 1L]

ik_ff <- unique(data[, .(gvkey, fyear, capx, ppegt)])
setorder(ik_ff, gvkey, fyear)
ik_ff[, ppegt_lag := shift(ppegt), by = gvkey]
ik_ff[, ik_ff := fifelse(
  !is.na(capx) & !is.na(ppegt_lag) & ppegt_lag > 0,
  capx / ppegt_lag,
  NA_real_
)]
ik_ff[, fyear := fyear + 1L]

### 9.8 June Formation Table ##################################################
june_form <- data[
  month == 6 &
    has_dec == TRUE &
    has_june == TRUE,
  .(
    permno,
    gvkey,
    ffyear = ffyear + 1L,
    exchcd,
    me_june = me_clean,
    bm = bm_ff,
    chi_sic = chi_sic
  )
]

june_form <- merge(
  june_form,
  op_ff[, .(gvkey, fyear, op_ff)],
  by.x = c("gvkey", "ffyear"),
  by.y = c("gvkey", "fyear"),
  all.x = TRUE
)

june_form <- merge(
  june_form,
  inv_ff[, .(gvkey, fyear, inv_ff)],
  by.x = c("gvkey", "ffyear"),
  by.y = c("gvkey", "fyear"),
  all.x = TRUE
)

june_form <- merge(
  june_form,
  ik_ff[, .(gvkey, fyear, ik_ff)],
  by.x = c("gvkey", "ffyear"),
  by.y = c("gvkey", "fyear"),
  all.x = TRUE
)

setnames(
  june_form,
  old = c("op_ff", "inv_ff", "ik_ff"),
  new = c("op", "inv", "ik")
)

setorder(june_form, permno, ffyear)
june_form <- unique(june_form, by = c("permno", "ffyear"))

### 9.9 NYSE Breakpoints for BM ###############################################
signals_nyse <- june_form[exchcd == 1]

bm_bp <- signals_nyse[
  !is.na(bm),
  {
    b <- ff_30_70_breaks(bm)
    .(bm30 = b["low"], bm70 = b["high"])
  },
  by = ffyear
]

june_form <- merge(
  june_form,
  bm_bp,
  by = "ffyear",
  all.x = TRUE
)

### 9.10 Assign HML Legs ######################################################
june_form[, value_leg := fcase(
  !is.na(bm) & bm <= bm30, "L",
  !is.na(bm) & bm > bm70, "H",
  default = NA_character_
)]

### 9.11 Return Panel: HML, CMA, RMW, and Baseline SIC SMQ ####################
return_panel <- merge(
  ff5[, .(mdate, HML, CMA, RMW)],
  SMQ[, .(mdate, SMQ_sic_VW)],
  by = "mdate",
  all = FALSE
)

return_panel <- return_panel[
  mdate >= START_DATE &
    mdate <= END_DATE
]

return_panel[, period := fifelse(
  mdate <= as.Date("2002-12-01"),
  "Early",
  "Late"
)]

### 9.12 Subperiod Summary Statistics #########################################
subperiod_stats <- rbindlist(
  lapply(names(subperiods), function(nm) {
    sp <- subperiods[[nm]]
    dt <- return_panel[
      mdate >= sp$start &
        mdate <= sp$end
    ]

    hml_test <- nw_mean_test(dt$HML, NW_LAGS)
    cma_test <- nw_mean_test(dt$CMA, NW_LAGS)
    rmw_test <- nw_mean_test(dt$RMW, NW_LAGS)
    smq_test <- nw_mean_test(dt$SMQ_sic_VW, NW_LAGS)

    data.table(
      period    = nm,
      HML_mean  = mean(dt$HML, na.rm = TRUE),
      HML_sd    = sd(dt$HML, na.rm = TRUE),
      HML_sr    = mean(dt$HML, na.rm = TRUE) / sd(dt$HML, na.rm = TRUE) * sqrt(12),
      HML_t     = hml_test$tstat,
      HML_p     = hml_test$pval,
      HML_s     = hml_test$stars,
      CMA_mean  = mean(dt$CMA, na.rm = TRUE),
      CMA_sd    = sd(dt$CMA, na.rm = TRUE),
      CMA_sr    = mean(dt$CMA, na.rm = TRUE) / sd(dt$CMA, na.rm = TRUE) * sqrt(12),
      CMA_t     = cma_test$tstat,
      CMA_p     = cma_test$pval,
      CMA_s     = cma_test$stars,
      RMW_mean  = mean(dt$RMW, na.rm = TRUE),
      RMW_sd    = sd(dt$RMW, na.rm = TRUE),
      RMW_sr    = mean(dt$RMW, na.rm = TRUE) / sd(dt$RMW, na.rm = TRUE) * sqrt(12),
      RMW_t     = rmw_test$tstat,
      RMW_p     = rmw_test$pval,
      RMW_s     = rmw_test$stars,
      SMQ_mean  = mean(dt$SMQ_sic_VW, na.rm = TRUE),
      SMQ_sd    = sd(dt$SMQ_sic_VW, na.rm = TRUE),
      SMQ_sr    = mean(dt$SMQ_sic_VW, na.rm = TRUE) / sd(dt$SMQ_sic_VW, na.rm = TRUE) * sqrt(12),
      SMQ_t     = smq_test$tstat,
      SMQ_p     = smq_test$pval,
      SMQ_s     = smq_test$stars
    )
  }),
  fill = TRUE
)

### 9.13 Subperiod Correlations ###############################################
subperiod_corrs <- list()

for (nm in names(subperiods)) {
  sp <- subperiods[[nm]]

  corr <- cor(
    return_panel[
      mdate >= sp$start &
        mdate <= sp$end,
      .(HML, CMA, RMW, SMQ_sic_VW)
    ],
    use = "pairwise.complete.obs"
  )

  subperiod_corrs[[nm]] <- corr
}

### 9.14 Conditional HML Regressions by Subperiod #############################
cond_regs <- list()

for (nm in names(subperiods)) {
  sp <- subperiods[[nm]]

  dt_p <- return_panel[
    mdate >= sp$start &
      mdate <= sp$end,
    .(HML, CMA, RMW, SMQ_sic_VW)
  ]

  reg1 <- run_ts_reg_nw(HML ~ CMA, dt_p, NW_LAGS)
  reg2 <- run_ts_reg_nw(HML ~ RMW, dt_p, NW_LAGS)
  reg3 <- run_ts_reg_nw(HML ~ SMQ_sic_VW, dt_p, NW_LAGS)
  reg4 <- run_ts_reg_nw(HML ~ CMA + RMW, dt_p, NW_LAGS)
  reg5 <- run_ts_reg_nw(HML ~ CMA + SMQ_sic_VW, dt_p, NW_LAGS)
  reg6 <- run_ts_reg_nw(HML ~ RMW + SMQ_sic_VW, dt_p, NW_LAGS)
  reg7 <- run_ts_reg_nw(HML ~ CMA + RMW + SMQ_sic_VW, dt_p, NW_LAGS)

  if (!is.null(reg1)) reg1[, `:=`(period = nm, spec = "HML ~ CMA")]
  if (!is.null(reg2)) reg2[, `:=`(period = nm, spec = "HML ~ RMW")]
  if (!is.null(reg3)) reg3[, `:=`(period = nm, spec = "HML ~ SMQ_sic_VW")]
  if (!is.null(reg4)) reg4[, `:=`(period = nm, spec = "HML ~ CMA + RMW")]
  if (!is.null(reg5)) reg5[, `:=`(period = nm, spec = "HML ~ CMA + SMQ_sic_VW")]
  if (!is.null(reg6)) reg6[, `:=`(period = nm, spec = "HML ~ RMW + SMQ_sic_VW")]
  if (!is.null(reg7)) reg7[, `:=`(period = nm, spec = "HML ~ CMA + RMW + SMQ_sic_VW")]

  cond_regs[[paste0(nm, "_1")]] <- reg1
  cond_regs[[paste0(nm, "_2")]] <- reg2
  cond_regs[[paste0(nm, "_3")]] <- reg3
  cond_regs[[paste0(nm, "_4")]] <- reg4
  cond_regs[[paste0(nm, "_5")]] <- reg5
  cond_regs[[paste0(nm, "_6")]] <- reg6
  cond_regs[[paste0(nm, "_7")]] <- reg7
}

conditional_table <- rbindlist(
  cond_regs[!sapply(cond_regs, is.null)],
  fill = TRUE
)

setcolorder(
  conditional_table,
  c("period", "spec", "term", "estimate", "tstat", "pval", "stars", "adj_r2")
)

conditional_alpha_table <- conditional_table[
  term == "(Intercept)",
  .(
    period,
    spec,
    alpha    = estimate,
    t_alpha  = tstat,
    p_alpha  = pval,
    alpha_s  = stars,
    adj_r2
  )
]

### Additional SMQ Spanning Tests by Subperiod ###############################
smq_spanning_specs <- list(
  "SMQ ~ CMA"             = SMQ_sic_VW ~ CMA,
  "SMQ ~ HML"             = SMQ_sic_VW ~ HML,
  "SMQ ~ RMW"             = SMQ_sic_VW ~ RMW,
  "SMQ ~ HML + CMA"       = SMQ_sic_VW ~ HML + CMA,
  "SMQ ~ HML + RMW"       = SMQ_sic_VW ~ HML + RMW,
  "SMQ ~ CMA + RMW"       = SMQ_sic_VW ~ CMA + RMW,
  "SMQ ~ HML + CMA + RMW" = SMQ_sic_VW ~ HML + CMA + RMW
)

smq_spanning_regs <- list()

for (nm in names(subperiods)) {
  sp <- subperiods[[nm]]

  dt_p <- return_panel[
    mdate >= sp$start &
      mdate <= sp$end,
    .(SMQ_sic_VW, HML, CMA, RMW)
  ]

  for (spec_name in names(smq_spanning_specs)) {
    reg <- run_ts_reg_nw(smq_spanning_specs[[spec_name]], dt_p, NW_LAGS)

    if (!is.null(reg)) {
      reg[, `:=`(period = nm, spec = spec_name)]
      smq_spanning_regs[[paste0(nm, "_", spec_name)]] <- reg
    }
  }
}

smq_spanning_table <- rbindlist(
  smq_spanning_regs[!sapply(smq_spanning_regs, is.null)],
  fill = TRUE
)

setcolorder(
  smq_spanning_table,
  c("period", "spec", "term", "estimate", "tstat", "pval", "stars", "adj_r2")
)

smq_spanning_alpha_table <- smq_spanning_table[
  term == "(Intercept)",
  .(
    period,
    spec,
    alpha    = estimate,
    t_alpha  = tstat,
    p_alpha  = pval,
    alpha_s  = stars,
    adj_r2
  )
]

### 9.15 Economic Composition of HML by Subperiod #############################
hml_chars_annual <- june_form[
  !is.na(value_leg) &
    value_leg %in% c("H", "L") &
    !is.na(me_june) &
    me_june > 0,
  .(
    inv_leg = weighted.mean(inv, me_june, na.rm = TRUE),
    ik_leg  = weighted.mean(ik, me_june, na.rm = TRUE),
    n_firms = uniqueN(permno)
  ),
  by = .(ffyear, value_leg)
]

hml_chars_panel <- merge(
  hml_chars_annual[value_leg == "H"],
  hml_chars_annual[value_leg == "L"],
  by = "ffyear",
  suffixes = c("_H", "_L"),
  all = FALSE
)

hml_chars_panel[, `:=`(
  spread_inv = inv_leg_H - inv_leg_L,
  spread_ik  = ik_leg_H - ik_leg_L
)]

hml_chars_panel[, ffyear := as.integer(ffyear)]
hml_chars_panel[, mdate := as.Date(sprintf("%04d-07-01", ffyear))]

hml_chars_panel <- hml_chars_panel[
  mdate >= START_DATE &
    mdate <= END_DATE
]

hml_chars_panel[, period := fifelse(
  mdate <= as.Date("2002-12-01"),
  "Early",
  "Late"
)]

leg_characteristics <- hml_chars_panel[
  ,
  .(
    H_inv   = mean(inv_leg_H, na.rm = TRUE),
    L_inv   = mean(inv_leg_L, na.rm = TRUE),
    HmL_inv = mean(spread_inv, na.rm = TRUE),
    H_ik    = mean(ik_leg_H, na.rm = TRUE),
    L_ik    = mean(ik_leg_L, na.rm = TRUE),
    HmL_ik  = mean(spread_ik, na.rm = TRUE)
  ),
  by = period
]

### 9.16 Mean Tests for HML Leg Investment Spreads ############################
spread_tests <- rbindlist(
  lapply(names(subperiods), function(nm) {
    sp <- subperiods[[nm]]
    dt <- hml_chars_panel[
      mdate >= sp$start &
        mdate <= sp$end
    ]

    inv_test <- nw_mean_test(dt$spread_inv, NW_LAGS)
    ik_test <- nw_mean_test(dt$spread_ik, NW_LAGS)

    data.table(
      period   = nm,
      variable = c("spread_inv", "spread_ik"),
      estimate = c(inv_test$estimate, ik_test$estimate),
      tstat    = c(inv_test$tstat, ik_test$tstat),
      pval     = c(inv_test$pval, ik_test$pval),
      stars    = c(inv_test$stars, ik_test$stars)
    )
  }),
  fill = TRUE
)

### 9.17 Compact Summary ######################################################
summary_compact <- data.table(
  Statistic = c(
    "HML mean early",
    "HML mean late",
    "CMA mean early",
    "CMA mean late",
    "RMW mean early",
    "RMW mean late",
    "SMQ_sic_VW mean early",
    "SMQ_sic_VW mean late",
    "H-L inv early",
    "H-L inv late",
    "H-L ik early",
    "H-L ik late"
  ),
  Value = c(
    subperiod_stats[period == "Early", HML_mean],
    subperiod_stats[period == "Late", HML_mean],
    subperiod_stats[period == "Early", CMA_mean],
    subperiod_stats[period == "Late", CMA_mean],
    subperiod_stats[period == "Early", RMW_mean],
    subperiod_stats[period == "Late", RMW_mean],
    subperiod_stats[period == "Early", SMQ_mean],
    subperiod_stats[period == "Late", SMQ_mean],
    leg_characteristics[period == "Early", HmL_inv],
    leg_characteristics[period == "Late", HmL_inv],
    leg_characteristics[period == "Early", HmL_ik],
    leg_characteristics[period == "Late", HmL_ik]
  ),
  tstat = c(
    subperiod_stats[period == "Early", HML_t],
    subperiod_stats[period == "Late", HML_t],
    subperiod_stats[period == "Early", CMA_t],
    subperiod_stats[period == "Late", CMA_t],
    subperiod_stats[period == "Early", RMW_t],
    subperiod_stats[period == "Late", RMW_t],
    subperiod_stats[period == "Early", SMQ_t],
    subperiod_stats[period == "Late", SMQ_t],
    spread_tests[period == "Early" & variable == "spread_inv", tstat],
    spread_tests[period == "Late" & variable == "spread_inv", tstat],
    spread_tests[period == "Early" & variable == "spread_ik", tstat],
    spread_tests[period == "Late" & variable == "spread_ik", tstat]
  )
)

### 9.18 Rolling Means and Rolling Sharpe Ratios ##############################
plot_panel <- copy(return_panel)

plot_panel[, HML_roll := roll_mean_safe(HML, ROLL_WIN)]
plot_panel[, CMA_roll := roll_mean_safe(CMA, ROLL_WIN)]
plot_panel[, RMW_roll := roll_mean_safe(RMW, ROLL_WIN)]
plot_panel[, SMQ_roll := roll_mean_safe(SMQ_sic_VW, ROLL_WIN)]

plot_panel[, HML_sr_roll := roll_sharpe_safe(HML, ROLL_WIN)]
plot_panel[, CMA_sr_roll := roll_sharpe_safe(CMA, ROLL_WIN)]
plot_panel[, RMW_sr_roll := roll_sharpe_safe(RMW, ROLL_WIN)]
plot_panel[, SMQ_sr_roll := roll_sharpe_safe(SMQ_sic_VW, ROLL_WIN)]

### 9.19 Diagnostics ##########################################################
cat("\n================ DIAGNOSTICS ================\n")
cat("Rows in june_form:", nrow(june_form), "\n")
cat("Rows in hml_chars_annual:", nrow(hml_chars_annual), "\n")
cat("Rows in hml_chars_panel:", nrow(hml_chars_panel), "\n")
cat("Rows in return_panel:", nrow(return_panel), "\n")

cat("\nValue-leg counts in june_form:\n")
print(june_form[, .N, by = value_leg][order(value_leg)])

cat("\nCharacteristic availability in hml_chars_panel:\n")
print(hml_chars_panel[, .(
  n_inv = sum(!is.na(spread_inv)),
  n_ik  = sum(!is.na(spread_ik))
)])

### 9.20 Plot Rolling Means ###################################################
par(mfrow = c(1, 1))
plot(
  plot_panel$mdate,
  plot_panel$HML_roll,
  type = "l",
  lwd = 2,
  main = "60-Month Rolling Means: HML, CMA, RMW, and Baseline SMQ_sic_VW",
  xlab = "",
  ylab = "Rolling Mean Return"
)
lines(plot_panel$mdate, plot_panel$CMA_roll, lwd = 2, lty = 2)
lines(plot_panel$mdate, plot_panel$RMW_roll, lwd = 2, lty = 4)
lines(plot_panel$mdate, plot_panel$SMQ_roll, lwd = 2, lty = 3)
abline(v = as.Date("2003-01-01"), lty = 5)
legend(
  "topright",
  legend = c("HML", "CMA", "RMW", "SMQ_sic_VW", "Early/Late split"),
  lty = c(1, 2, 4, 3, 5),
  lwd = c(2, 2, 2, 2, 1),
  bty = "n"
)

### 9.21 Plot Rolling Sharpe Ratios ###########################################
par(mfrow = c(1, 1))
plot(
  plot_panel$mdate,
  plot_panel$HML_sr_roll,
  type = "l",
  lwd = 2,
  main = "60-Month Rolling Sharpe Ratios: HML, CMA, RMW, and Baseline SMQ_sic_VW",
  xlab = "",
  ylab = "Rolling Sharpe Ratio"
)
lines(plot_panel$mdate, plot_panel$CMA_sr_roll, lwd = 2, lty = 2)
lines(plot_panel$mdate, plot_panel$RMW_sr_roll, lwd = 2, lty = 4)
lines(plot_panel$mdate, plot_panel$SMQ_sr_roll, lwd = 2, lty = 3)
abline(v = as.Date("2003-01-01"), lty = 5)
legend(
  "topright",
  legend = c("HML", "CMA", "RMW", "SMQ_sic_VW", "Early/Late split"),
  lty = c(1, 2, 4, 3, 5),
  lwd = c(2, 2, 2, 2, 1),
  bty = "n"
)

### 9.22 Save #################################################################
saveRDS(return_panel, "code9_return_panel_subperiods_sic_baseline.rds")
saveRDS(subperiod_stats, "code9_subperiod_stats_investment_q_sic_baseline.rds")
saveRDS(subperiod_corrs, "code9_subperiod_correlations_investment_q_sic_baseline.rds")
saveRDS(conditional_table, "code9_conditional_hml_regs_investment_q_sic_baseline.rds")
saveRDS(conditional_alpha_table, "code9_conditional_hml_alphas_investment_q_sic_baseline.rds")
saveRDS(hml_chars_panel, "code9_hml_leg_chars_panel_investment_q_sic_baseline.rds")
saveRDS(leg_characteristics, "code9_hml_leg_characteristics_by_period_sic_baseline.rds")
saveRDS(spread_tests, "code9_hml_investment_spread_tests_sic_baseline.rds")
saveRDS(summary_compact, "code9_summary_compact_investment_q_sic_baseline.rds")
saveRDS(plot_panel, "code9_plot_panel_with_rolling_mean_and_sharpe_sic_baseline.rds")
saveRDS(smq_spanning_table, "code9_smq_spanning_regs_by_subperiod_sic_baseline.rds")
saveRDS(smq_spanning_alpha_table, "code9_smq_spanning_alphas_by_subperiod_sic_baseline.rds")

### 9.23 Print ################################################################
cat("\n================ SUBPERIOD SUMMARY STATISTICS ================\n")
print(round_dt(subperiod_stats, 4))

cat("\n================ SUBPERIOD CORRELATIONS: EARLY ================\n")
print(round(subperiod_corrs$Early, 4))

cat("\n================ SUBPERIOD CORRELATIONS: LATE ================\n")
print(round(subperiod_corrs$Late, 4))

cat("\n================ CONDITIONAL HML REGRESSIONS ================\n")
print(round_dt(conditional_table, 4))

cat("\n================ CONDITIONAL HML ALPHAS ONLY ================\n")
print(round_dt(conditional_alpha_table, 4))

cat("\n================ SMQ SPANNING REGRESSIONS BY SUBPERIOD ================\n")
print(round_dt(smq_spanning_table, 4))

cat("\n================ SMQ SPANNING ALPHAS BY SUBPERIOD ================\n")
print(round_dt(smq_spanning_alpha_table, 4))

cat("\n================ HML LEG CHARACTERISTICS BY SUBPERIOD ================\n")
print(round_dt(leg_characteristics, 4))

cat("\n================ HML INVESTMENT-SPREAD TESTS ================\n")
print(round_dt(spread_tests, 4))

cat("\n================ COMPACT SUMMARY ================\n")
print(round_dt(summary_compact, 4))

cat("\nRolling means are computed as backward-looking 60-month averages.\n")
cat("Rolling Sharpe ratios are computed as backward-looking 60-month mean/sd * sqrt(12).\n")




### 9.14 HML and CMA Summary and Spanning Tests ##############################

### Separate HML-CMA panel
### Do not merge with SMQ, since this would truncate the sample
hml_cma_panel <- ff5[
  mdate >= as.Date("1963-07-01") &
    mdate <= as.Date("2024-12-01"),
  .(
    mdate,
    HML,
    CMA
  )
]

setorder(hml_cma_panel, mdate)

### Remove months missing either HML or CMA
### This ensures that summary statistics and regressions use the same sample
hml_cma_panel <- hml_cma_panel[
  complete.cases(HML, CMA)
]

### Define full period and subperiods
hml_cma_periods <- list(
  Full = list(
    start = as.Date("1963-07-01"),
    end   = as.Date("2024-12-01")
  ),
  Early = list(
    start = as.Date("1963-07-01"),
    end   = as.Date("1993-12-01")
  ),
  Late = list(
    start = as.Date("1994-01-01"),
    end   = as.Date("2024-12-01")
  )
)

### Estimate summary statistics and spanning regressions
hml_cma_results <- rbindlist(
  lapply(names(hml_cma_periods), function(nm) {
    sp <- hml_cma_periods[[nm]]

    dt_p <- hml_cma_panel[
      mdate >= sp$start &
        mdate <= sp$end
    ]

    if (nrow(dt_p) < 24) {
      stop(
        paste(
          "Fewer than 24 observations available for",
          nm,
          "period."
        )
      )
    }

    ### Newey-West mean-return tests
    hml_mean_test <- nw_mean_test(
      dt_p$HML,
      NW_LAGS
    )

    cma_mean_test <- nw_mean_test(
      dt_p$CMA,
      NW_LAGS
    )

    ### Pairwise spanning regressions
    ###
    ### HML_t = alpha + beta CMA_t + error_t
    hml_on_cma <- run_ts_reg_nw(
      HML ~ CMA,
      dt_p,
      NW_LAGS
    )

    ### CMA_t = alpha + beta HML_t + error_t
    cma_on_hml <- run_ts_reg_nw(
      CMA ~ HML,
      dt_p,
      NW_LAGS
    )

    ### Extract intercepts
    hml_alpha <- hml_on_cma[
      term == "(Intercept)"
    ]

    cma_alpha <- cma_on_hml[
      term == "(Intercept)"
    ]

    ### Extract factor loadings
    hml_beta <- hml_on_cma[
      term == "CMA"
    ]

    cma_beta <- cma_on_hml[
      term == "HML"
    ]

    data.table(
      period = nm,
      start_date = min(dt_p$mdate),
      end_date = max(dt_p$mdate),
      factor = c(
        "HML",
        "CMA"
      ),
      spanning_factor = c(
        "CMA",
        "HML"
      ),
      mean = c(
        mean(dt_p$HML),
        mean(dt_p$CMA)
      ),
      t_mean = c(
        hml_mean_test$tstat,
        cma_mean_test$tstat
      ),
      p_mean = c(
        hml_mean_test$pval,
        cma_mean_test$pval
      ),
      mean_stars = c(
        hml_mean_test$stars,
        cma_mean_test$stars
      ),
      sd = c(
        sd(dt_p$HML),
        sd(dt_p$CMA)
      ),
      sharpe = c(
        mean(dt_p$HML) /
          sd(dt_p$HML) * sqrt(12),
        mean(dt_p$CMA) /
          sd(dt_p$CMA) * sqrt(12)
      ),
      alpha = c(
        hml_alpha$estimate,
        cma_alpha$estimate
      ),
      t_alpha = c(
        hml_alpha$tstat,
        cma_alpha$tstat
      ),
      p_alpha = c(
        hml_alpha$pval,
        cma_alpha$pval
      ),
      alpha_stars = c(
        hml_alpha$stars,
        cma_alpha$stars
      ),
      beta = c(
        hml_beta$estimate,
        cma_beta$estimate
      ),
      t_beta = c(
        hml_beta$tstat,
        cma_beta$tstat
      ),
      p_beta = c(
        hml_beta$pval,
        cma_beta$pval
      ),
      beta_stars = c(
        hml_beta$stars,
        cma_beta$stars
      ),
      adj_r2 = c(
        unique(hml_on_cma$adj_r2)[1],
        unique(cma_on_hml$adj_r2)[1]
      ),
      n_months = nrow(dt_p)
    )
  }),
  fill = TRUE
)

### Arrange rows
hml_cma_results[
  ,
  period_order := match(
    period,
    c("Full", "Early", "Late")
  )
]

hml_cma_results[
  ,
  factor_order := match(
    factor,
    c("HML", "CMA")
  )
]

setorder(
  hml_cma_results,
  period_order,
  factor_order
)

hml_cma_results[
  ,
  c(
    "period_order",
    "factor_order"
  ) := NULL
]

### Presentation table
hml_cma_display <- hml_cma_results[
  ,
  .(
    Period = period,
    `Sample period` = paste0(
      format(start_date, "%Y-%m"),
      " to ",
      format(end_date, "%Y-%m")
    ),
    Factor = factor,
    `Mean excess return (%)` = paste0(
      sprintf("%.3f", 100 * mean),
      mean_stars
    ),
    `t-statistic mean` = sprintf(
      "%.2f",
      t_mean
    ),
    `Standard deviation (%)` = sprintf(
      "%.3f",
      100 * sd
    ),
    `Annualized Sharpe ratio` = sprintf(
      "%.3f",
      sharpe
    ),
    `Spanned by` = spanning_factor,
    `Alpha (%)` = paste0(
      sprintf("%.3f", 100 * alpha),
      alpha_stars
    ),
    `t-statistic alpha` = sprintf(
      "%.2f",
      t_alpha
    ),
    Loading = paste0(
      sprintf("%.3f", beta),
      beta_stars
    ),
    `t-statistic loading` = sprintf(
      "%.2f",
      t_beta
    ),
    `Adjusted R2` = sprintf(
      "%.3f",
      adj_r2
    ),
    Observations = n_months
  )
]

### Diagnostics
cat("\n================ HML-CMA SAMPLE DIAGNOSTICS ================\n")

cat(
  "First available month:",
  format(min(hml_cma_panel$mdate), "%Y-%m"),
  "\n"
)

cat(
  "Last available month:",
  format(max(hml_cma_panel$mdate), "%Y-%m"),
  "\n"
)

cat(
  "Total complete observations:",
  nrow(hml_cma_panel),
  "\n"
)

### Print final table
cat("\n================ HML AND CMA SPANNING RESULTS ================\n")

print(hml_cma_display)
