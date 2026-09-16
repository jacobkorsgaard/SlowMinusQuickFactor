## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026

### 2. Empirical Variable Construction ########################################
### 2.1 Setup ##################################################################
cat("\014")
rm(list = ls())
graphics.off()

library(data.table)
library(lubridate)

source(file.path("R", "00_config.R"))

### 2.2 Load Data ##############################################################
data <- readRDS("crsp_compustat_monthly_panel_1950_2024.rds")
setDT(data)

## Ensure SIC variables are numeric and available throughout
data[, sic := as.integer(sic)]
data[, sic2 := fifelse(!is.na(sic), sic %/% 100L, NA_integer_)]

## Basic checks
names(data)
range(data$mdate, na.rm = TRUE)

### 2.3 Calendar Definitions ###################################################
data[, year := year(mdate)]
data[, month := month(mdate)]
data[, ffyear := fifelse(month >= 7, year, year - 1)]

### 2.4 Sample Eligibility Requirements ########################################
### 2.4.1 Compustat 2-Year Requirement #########################################
## At June formation, require that fiscal years ffyear - 1 and ffyear - 2 exist.
## Given your ffyear convention, June of calendar year t has ffyear = t - 1.
## So this condition mirrors your earlier logic and is kept intact.
comp_ok <-
  data[
    month == 6,
    .(
      comp_ok =
        any(fyear == ffyear - 1) &
          any(fyear == ffyear - 2)
    ),
    by = .(gvkey, ffyear)
  ]

data <-
  merge(
    data,
    comp_ok,
    by = c("gvkey", "ffyear"),
    all.x = TRUE
  )

### 2.4.2 CRSP Presence Requirement ############################################
## Require presence in December of t - 1 and June of t for June t formation.
dec_presence <-
  data[
    month == 12,
    .(has_dec = TRUE),
    by = .(permno, ffyear = year + 1L)
  ]

june_presence <-
  data[
    month == 6,
    .(has_june = TRUE),
    by = .(permno, ffyear)
  ]

data <-
  merge(
    data,
    dec_presence,
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

data <-
  merge(
    data,
    june_presence,
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

### 2.5 Returns and Prices #####################################################
## Use absolute price because CRSP stores negative prices as bid/ask conventions.
data[, prc := abs(prc)]

## Delisting-adjusted return
data[, retadj :=
  fifelse(
    !is.na(ret) & !is.na(dlret),
    (1 + ret) * (1 + dlret) - 1,
    fifelse(
      !is.na(ret),
      ret,
      fifelse(!is.na(dlret), dlret, NA_real_)
    )
  )]

### 2.6 Market Equity ##########################################################
## CRSP market equity
data[, mve :=
  fifelse(
    !is.na(prc) & !is.na(shrout) & shrout > 0,
    prc * shrout,
    NA_real_
  )]

## Keep a clean version in original CRSP units
data[, me_clean := mve]

## Convert to millions to align with Compustat scale
data[, me_m := me_clean / 1000]

## December market equity aligned to holding ffyear
data[, me_dec := fifelse(month == 12, me_m, NA_real_)]

## Carry December t - 1 ME into ffyear = t
data[, me_dec_cf := {
  x <- me_dec
  if (length(x[!is.na(x)]) == 1) x[!is.na(x)][1] else NA_real_
}, by = .(permno, ffyear)]

### 2.7 Book Equity ############################################################
## Preferred stock hierarchy: redemption > liquidation > carrying value
data[, ps :=
  fifelse(
    !is.na(pstkrv),
    pstkrv,
    fifelse(
      !is.na(pstkl),
      pstkl,
      fifelse(!is.na(pstk), pstk, 0)
    )
  )]

## Fama-French book equity
data[, be :=
  fifelse(
    !is.na(seq),
    seq + fifelse(!is.na(txditc), txditc, 0) - ps,
    NA_real_
  )]

## Keep positive BE only
data[, be_clean :=
  fifelse(
    !is.na(be) & be > 0,
    be,
    NA_real_
  )]

### 2.8 Book-to-Market #########################################################
data[, bm_ff :=
  fifelse(
    !is.na(be_clean) &
      !is.na(me_dec_cf) &
      me_dec_cf > 0,
    be_clean / me_dec_cf,
    NA_real_
  )]

data[, bm_ff := fifelse(bm_ff > 0, bm_ff, NA_real_)]

### 2.9 Investment #############################################################
## Build annual asset panel
annual_at <- unique(data[, .(gvkey, fyear, at)], by = c("gvkey", "fyear"))
setorder(annual_at, gvkey, fyear)

annual_at[, at_lag := shift(at), by = gvkey]
annual_at[, inv_annual :=
  fifelse(
    !is.na(at) & !is.na(at_lag) & at_lag > 0,
    (at - at_lag) / at_lag,
    NA_real_
  )]

data <-
  merge(
    data,
    annual_at[, .(gvkey, fyear, inv_annual)],
    by = c("gvkey", "fyear"),
    all.x = TRUE
  )

### 2.10 Operating Profitability ###############################################
## Keep your strict version: require all major components to be observed.
data[, op :=
  fifelse(
    !is.na(revt) &
      !is.na(cogs) &
      !is.na(xsga) &
      !is.na(xint) &
      !is.na(be_clean),
    (revt - cogs - xsga - xint) / be_clean,
    NA_real_
  )]

data[, op_clean := fifelse(!is.na(op), op, NA_real_)]

### 2.11 Physical Investment (I/K) #############################################
annual_ik <-
  data[
    !is.na(capx) & !is.na(ppegt) & ppegt > 0,
    .(
      capx = first(capx),
      ppegt = first(ppegt)
    ),
    by = .(gvkey, fyear)
  ]

setorder(annual_ik, gvkey, fyear)

annual_ik[, ppegt_lag := shift(ppegt), by = gvkey]
annual_ik[, ik_annual :=
  fifelse(
    !is.na(capx) & !is.na(ppegt_lag) & ppegt_lag > 0,
    capx / ppegt_lag,
    NA_real_
  )]

data <-
  merge(
    data,
    annual_ik[, .(gvkey, fyear, ppegt_lag, ik_annual)],
    by = c("gvkey", "fyear"),
    all.x = TRUE
  )

### 2.12 Tobin's Q #############################################################
## Compustat market equity
data[, mve_cs :=
  fifelse(
    !is.na(prcc_f) & !is.na(csho) & csho > 0,
    prcc_f * csho,
    NA_real_
  )]

## Enterprise-value style numerator
data[, firm_value_cs :=
  fifelse(
    !is.na(mve_cs),
    mve_cs +
      fifelse(!is.na(dltt), dltt, 0) +
      fifelse(!is.na(dlc), dlc, 0) -
      fifelse(!is.na(act), act, 0),
    NA_real_
  )]

## Tobin's q
data[, tobins_q :=
  fifelse(
    !is.na(firm_value_cs) & !is.na(ppegt) & ppegt > 5,
    firm_value_cs / ppegt,
    NA_real_
  )]

### 2.13 June Formation Table ##################################################
## Build one June record per permno holding year.
## June of calendar year t is mapped into holding ffyear = t.
june_form <-
  data[
    month == 6 &
      has_dec == TRUE &
      has_june == TRUE,
    .(
      gvkey = gvkey[1],
      ffyear = ffyear[1] + 1L,
      bm_june = bm_ff[1],
      op_june = op_clean[1],
      inv_june = inv_annual[1],
      ik_june = ik_annual[1],
      q_june = tobins_q[1],
      me_june = me_clean[1],
      ppe_june = ppegt_lag[1]
    ),
    by = .(permno, ffyear)
  ]

## Reattach June signals to the monthly panel by permno and holding ffyear
data <-
  merge(
    data,
    june_form[, .(
      permno,
      ffyear,
      bm_cf = bm_june,
      op_cf = op_june,
      inv_cf = inv_june,
      ik_cf = ik_june,
      q_cf = q_june,
      me_ff = me_june,
      ppe_cf = ppe_june
    )],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

### 2.14 June Timing Diagnostics ###############################################
diag <-
  data[
    month %in% c(6, 7) &
      has_dec == TRUE &
      has_june == TRUE,
    .(
      n_rows = .N,
      n_bm_cf = sum(!is.na(bm_cf)),
      n_op_cf = sum(!is.na(op_cf)),
      n_inv_cf = sum(!is.na(inv_cf)),
      n_ik_cf = sum(!is.na(ik_cf)),
      n_q_cf = sum(!is.na(q_cf))
    ),
    by = month
  ]
print(diag)

check_july <-
  data[
    month %in% c(6, 7) &
      has_dec == TRUE &
      has_june == TRUE,
    .(
      june_bm = bm_ff[month == 6][1],
      july_bm_cf = bm_cf[month == 7][1]
    ),
    by = .(gvkey, ffyear)
  ]

print(check_july[!is.na(june_bm) & is.na(july_bm_cf)][1:20])

cat(
  "Share of gvkey-ffyear with June bm_ff but missing July bm_cf:",
  mean(!is.na(check_july$june_bm) & is.na(check_july$july_bm_cf), na.rm = TRUE),
  "\n"
)

### 2.15 Adjustment Cost Controls ##############################################
CHI_CTRL <- list(
  window_firm = 21L,
  window_group = 21L,
  min_obs_firm = 5L,
  min_obs_group = 5L,
  burn_in = 21L
)

### 2.16 Build Annual Estimation Panel #########################################
## Estimate q-investment relations using July snapshots of the carried-forward
## June signals. This matches your earlier logic, but now keeps true sic/sic2.
winsor_1pct <- function(x) {
  qs <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE, type = 7)
  x <- pmin(pmax(x, qs[1]), qs[2])
  x
}

annual <-
  data[
    month == 7 &
      !is.na(ik_cf) &
      !is.na(q_cf) &
      !is.na(sic) &
      !is.na(sic2) &
      !(sic >= 4900 & sic <= 4999) &
      !(sic >= 6000 & sic <= 6999) &
      !(sic >= 9000),
    .(
      invest = ik_cf[1],
      tobins_q = q_cf[1],
      sic = sic[1],
      sic2 = sic2[1]
    ),
    by = .(gvkey, fyear)
  ]

setorder(annual, gvkey, fyear)
annual <- unique(annual, by = c("gvkey", "fyear"))
setDT(annual)

annual[, invest := winsor_1pct(invest)]
annual[, tobins_q := winsor_1pct(tobins_q)]

annual[, q_lag := shift(tobins_q), by = gvkey]
annual <- annual[!is.na(q_lag)]

GLOBAL_START_YEAR <- min(annual$fyear)
BURN_IN <- CHI_CTRL$burn_in

### 2.17 Firm-Level Adjustment Costs ###########################################
### 2.17.1 Rolling ############################################################
estimate_firm_adjcost_roll <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    sub <- dt[fyear %in% ((t - CHI_CTRL$window_firm):(t - 1))]
    if (nrow(sub) < CHI_CTRL$min_obs_firm) next

    fit <- try(lm(invest ~ q_lag, data = sub), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        gvkey = unique(dt$gvkey),
        fyear = t,
        gamma_firm = g,
        chi_firm = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(sub)
      )
  }

  rbindlist(out, fill = TRUE)
}

regdata_firm <- annual[, .(gvkey, fyear, invest, q_lag)]
setorder(regdata_firm, gvkey, fyear)

adjcost_firm_panel <-
  regdata_firm[, estimate_firm_adjcost_roll(.SD), by = gvkey]

### 2.17.2 Expanding ##########################################################
estimate_firm_adjcost_exp <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    sub <- dt[fyear < t]
    if (nrow(sub) < CHI_CTRL$min_obs_firm) next

    fit <- try(lm(invest ~ q_lag, data = sub), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        gvkey = unique(dt$gvkey),
        fyear = t,
        gamma_firm_exp = g,
        chi_firm_exp = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(sub)
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_firm_exp_panel <-
  regdata_firm[, estimate_firm_adjcost_exp(.SD), by = gvkey]

### 2.17.2 Expanding ##########################################################
estimate_firm_adjcost_exp_YFE <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    sub <- dt[fyear < t]
    if (nrow(sub) < CHI_CTRL$min_obs_firm) next

    fit <- try(lm(invest ~ q_lag + factor(fyear), data = sub), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        gvkey = unique(dt$gvkey),
        fyear = t,
        gamma_firm_exp = g,
        chi_firm_exp = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(sub)
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_firm_exp_panel <-
  regdata_firm[, estimate_firm_adjcost_exp_YFE(.SD), by = gvkey]

### 2.18 SIC-Level Adjustment Costs ############################################
## Pooled within sic with firm FE and year FE
regdata_sic <-
  annual[
    !is.na(sic),
    .(gvkey, sic, fyear, invest, q_lag)
  ]
setorder(regdata_sic, sic, fyear, gvkey)

### 2.18.1 Rolling ############################################################
estimate_sic_adjcost_roll <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    win <- dt[fyear %in% ((t - CHI_CTRL$window_group):(t - 1))]
    if (nrow(win) < CHI_CTRL$min_obs_group) next

    fit <- try(lm(invest ~ q_lag + factor(gvkey) + factor(fyear), data = win), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        sic = unique(dt$sic),
        fyear = t,
        gamma_sic = g,
        chi_sic = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(win),
        adj_r2 = summary(fit)$adj.r.squared,
        se_gamma = sqrt(diag(vcov(fit)))[["q_lag"]]
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_sic_panel <-
  regdata_sic[, estimate_sic_adjcost_roll(.SD), by = sic]

### 2.18.2 Expanding ##########################################################
estimate_sic_adjcost_exp <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    sub <- dt[fyear < t]
    if (nrow(sub) < CHI_CTRL$min_obs_group) next

    fit <- try(lm(invest ~ q_lag + factor(gvkey) + factor(fyear), data = sub), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        sic = unique(dt$sic),
        fyear = t,
        gamma_sic_exp = g,
        chi_sic_exp = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(sub),
        adj_r2 = summary(fit)$adj.r.squared,
        se_gamma = sqrt(diag(vcov(fit)))[["q_lag"]]
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_sic_panel_exp <-
  regdata_sic[, estimate_sic_adjcost_exp(.SD), by = sic]

### 2.19 SIC2-Level Adjustment Costs ###########################################
## Pooled within sic2 with firm FE and year FE
regdata_sic2 <-
  annual[
    !is.na(sic2),
    .(gvkey, sic2, fyear, invest, q_lag)
  ]
setorder(regdata_sic2, sic2, fyear, gvkey)

### 2.19.1 Rolling ############################################################
estimate_sic2_adjcost_roll <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    win <- dt[fyear %in% ((t - CHI_CTRL$window_group):(t - 1))]
    if (nrow(win) < CHI_CTRL$min_obs_group) next

    fit <- try(lm(invest ~ q_lag + factor(gvkey) + factor(fyear), data = win), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        sic2 = unique(dt$sic2),
        fyear = t,
        gamma_sic2 = g,
        chi_sic2 = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(win),
        adj_r2 = summary(fit)$adj.r.squared,
        se_gamma = sqrt(diag(vcov(fit)))[["q_lag"]]
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_sic2_panel <-
  regdata_sic2[, estimate_sic2_adjcost_roll(.SD), by = sic2]

### 2.19.2 Expanding ##########################################################
estimate_sic2_adjcost_exp <- function(dt) {
  years <- sort(unique(dt$fyear))
  out <- list()

  for (t in years) {
    if (t < GLOBAL_START_YEAR + BURN_IN) next

    sub <- dt[fyear < t]
    if (nrow(sub) < CHI_CTRL$min_obs_group) next

    fit <- try(lm(invest ~ q_lag + factor(gvkey) + factor(fyear), data = sub), silent = TRUE)
    if (inherits(fit, "try-error")) next

    g <- coef(fit)[["q_lag"]]

    out[[as.character(t)]] <-
      data.table(
        sic2 = unique(dt$sic2),
        fyear = t,
        gamma_sic2_exp = g,
        chi_sic2_exp = ifelse(!is.na(g) & g != 0, 1 / g, NA_real_),
        n_obs = nrow(sub),
        adj_r2 = summary(fit)$adj.r.squared,
        se_gamma = sqrt(diag(vcov(fit)))[["q_lag"]]
      )
  }

  rbindlist(out, fill = TRUE)
}

adjcost_sic2_panel_exp <-
  regdata_sic2[, estimate_sic2_adjcost_exp(.SD), by = sic2]

### 2.20 Diagnostics ###########################################################
summary(adjcost_firm_panel$chi_firm)
summary(adjcost_firm_exp_panel$chi_firm_exp)

summary(adjcost_sic_panel$chi_sic)
summary(adjcost_sic_panel_exp$chi_sic_exp)

summary(adjcost_sic2_panel$chi_sic2)
summary(adjcost_sic2_panel_exp$chi_sic2_exp)

### 2.21 Merge Adjustment Cost Estimates #######################################
## Firm-level rolling
adjcost_firm_chi <- adjcost_firm_panel[, .(gvkey, fyear, chi_firm)]
setkey(data, gvkey, fyear)
setkey(adjcost_firm_chi, gvkey, fyear)
data <- adjcost_firm_chi[data]

## Firm-level expanding
adjcost_firm_chi_exp <- adjcost_firm_exp_panel[, .(gvkey, fyear, chi_firm_exp)]
setkey(adjcost_firm_chi_exp, gvkey, fyear)
data <- adjcost_firm_chi_exp[data]

## SIC-level rolling
adjcost_sic_chi <- adjcost_sic_panel[, .(sic, fyear, chi_sic)]
setkey(data, sic, fyear)
setkey(adjcost_sic_chi, sic, fyear)
data <- adjcost_sic_chi[data]

## SIC-level expanding
adjcost_sic_chi_exp <- adjcost_sic_panel_exp[, .(sic, fyear, chi_sic_exp)]
setkey(adjcost_sic_chi_exp, sic, fyear)
data <- adjcost_sic_chi_exp[data]

## SIC2-level rolling
adjcost_sic2_chi <- adjcost_sic2_panel[, .(sic2, fyear, chi_sic2)]
setkey(data, sic2, fyear)
setkey(adjcost_sic2_chi, sic2, fyear)
data <- adjcost_sic2_chi[data]

## SIC2-level expanding
adjcost_sic2_chi_exp <- adjcost_sic2_panel_exp[, .(sic2, fyear, chi_sic2_exp)]
setkey(adjcost_sic2_chi_exp, sic2, fyear)
data <- adjcost_sic2_chi_exp[data]

### 2.22 Final Variables and Checks ############################################
data[, exch :=
  fifelse(
    exchcd == 1,
    "NYSE",
    fifelse(
      exchcd == 2,
      "AMEX",
      fifelse(exchcd == 3, "NASDAQ", NA_character_)
    )
  )]

names(data)

summary(data$chi_firm)
summary(data$chi_firm_exp)
summary(data$chi_sic)
summary(data$chi_sic_exp)
summary(data$chi_sic2)
summary(data$chi_sic2_exp)

### 2.23 Save ##################################################################
saveRDS(data, "crsp_compustat_with_adj_cost.rds")
saveRDS(adjcost_firm_panel, "adjcost_firm_panel.rds")
saveRDS(adjcost_firm_exp_panel, "adjcost_firm_exp_panel.rds")
saveRDS(adjcost_sic_panel, "adjcost_sic_panel.rds")
saveRDS(adjcost_sic_panel_exp, "adjcost_sic_panel_exp.rds")
saveRDS(adjcost_sic2_panel, "adjcost_sic2_panel.rds")
saveRDS(adjcost_sic2_panel_exp, "adjcost_sic2_panel_exp.rds")
