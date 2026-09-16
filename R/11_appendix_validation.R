## A Appendix #################################################################
## A.0 Setup ##################################################################
cat("\014")
rm(list = ls())
graphics.off()

## Libraries
library(data.table)
library(lubridate)
library(sandwich)
library(lmtest)

## Project paths
source(file.path("R", "00_config.R"))

## Paths
APPENDIX_DIR <- INTERIM_DATA_DIR

## Core Files
PANEL_FILE <- "crsp_compustat_with_adj_cost.rds"
FF5_FILE <- "ff5_factors_internal.rds"
SMQ_FILE <- "SMQ_factors.rds"

## Sample Windows
START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")
PORT_START_DATE <- as.Date("1983-07-01")
PORT_END_DATE <- as.Date("2024-12-01")

## Newey-West Settings
NW_LAGS <- 12

## Load Data
data <- readRDS(PANEL_FILE)
setDT(data)
data[, mdate := as.Date(mdate)]

ff5 <- readRDS(FF5_FILE)
setDT(ff5)
ff5[, mdate := as.Date(mdate)]

SMQ <- readRDS(SMQ_FILE)
setDT(SMQ)
SMQ[, mdate := as.Date(mdate)]

## Refresh core calendar fields
data[, year := year(mdate)]
data[, month := month(mdate)]
data[, ffyear := fifelse(month >= 7, year, year - 1)]

ff5[, mdate := floor_date(mdate, "month")]
SMQ[, mdate := floor_date(mdate, "month")]

data[, sic := as.integer(sic)]

if (!"sic2" %in% names(data)) {
  data[, sic2 := fifelse(!is.na(sic), sic %/% 100L, NA_integer_)]
}

## Restrict panel window where relevant
data <- data[mdate >= START_DATE & mdate <= END_DATE]

## Helpers ####################################################################
filter_window <- function(dt, start_date, end_date) {
  dt[mdate >= start_date & mdate <= end_date]
}

cum_safe <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  if (any(ok)) {
    out[ok] <- cumprod(1 + x[ok]) - 1
  }
  out
}

winsor_1pct <- function(x) {
  q <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE, type = 7)
  pmin(pmax(x, q[1]), q[2])
}

star_p <- function(p) {
  fifelse(
    p < 0.01, "***",
    fifelse(
      p < 0.05, "**",
      fifelse(p < 0.10, "*", "")
    )
  )
}

round_dt <- function(dt, digits = 4) {
  out <- copy(dt)
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[, (num_cols) := lapply(.SD, round, digits), .SDcols = num_cols]
  out
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
  fit_sum <- summary(fit)

  out <- data.table(
    term      = rownames(ct),
    estimate  = as.numeric(ct[, 1]),
    tstat     = as.numeric(ct[, 3]),
    pval      = as.numeric(ct[, 4]),
    stars     = star_p(as.numeric(ct[, 4]))
  )

  out[, r2 := fit_sum$r.squared]
  out[, adj_r2 := fit_sum$adj.r.squared]
  out
}

get_reg_value <- function(x, term_name, col_name) {
  if (is.null(x) || !is.data.table(x) || nrow(x) == 0) {
    return(NA_real_)
  }
  if (!("term" %in% names(x))) {
    return(NA_real_)
  }

  out <- x[term == term_name, get(col_name)]
  if (length(out) == 0) {
    return(NA_real_)
  }
  as.numeric(out[1])
}

get_reg_star <- function(x, term_name) {
  if (is.null(x) || !is.data.table(x) || nrow(x) == 0) {
    return(NA_character_)
  }
  if (!("term" %in% names(x)) || !("stars" %in% names(x))) {
    return(NA_character_)
  }

  out <- x[term == term_name, stars]
  if (length(out) == 0) {
    return(NA_character_)
  }
  as.character(out[1])
}

## Helper: parse French 5x5 files #############################################
read_french_5x5 <- function(path, char_prefix = "V") {
  lines <- readLines(path, warn = FALSE)

  start_idx <- which(grepl("^\\s*\\d{6}\\s*,", lines))[1]
  if (is.na(start_idx)) stop("Could not find data start in: ", path)

  after <- lines[(start_idx + 1):length(lines)]
  blank_rel <- which(trimws(after) == "")[1]
  end_idx <- if (!is.na(blank_rel)) (start_idx + blank_rel - 1) else length(lines)

  dt <- fread(text = paste(lines[start_idx:end_idx], collapse = "\n"), header = FALSE)
  setDT(dt)

  if (ncol(dt) < 2) stop("French data block has too few columns in: ", path)

  setnames(dt, 1, "yyyymm")
  dt[, yyyymm := as.integer(yyyymm)]
  dt <- dt[!is.na(yyyymm)]

  dt[, year := yyyymm %/% 100L]
  dt[, mon := yyyymm %% 100L]
  dt[, mdate := as.Date(sprintf("%04d-%02d-01", year, mon))]
  dt[, c("year", "mon") := NULL]

  ret_cols <- setdiff(names(dt), c("yyyymm", "mdate"))
  for (cc in ret_cols) dt[, (cc) := as.numeric(get(cc)) / 100]

  dt[, yyyymm := NULL]

  ret_cols2 <- setdiff(names(dt), "mdate")
  if (length(ret_cols2) < 25) stop("Expected at least 25 return columns in: ", path)

  new_names <- paste0(
    rep(paste0("S", 1:5), each = 5),
    rep(paste0(char_prefix, 1:5), times = 5)
  )
  setnames(dt, old = ret_cols2[1:25], new = new_names)

  setcolorder(dt, c("mdate", setdiff(names(dt), "mdate")))
  dt
}

## Helper: parse French 2x4x4 files ###########################################
read_french_2x4x4 <- function(path) {
  lines <- readLines(path, warn = FALSE)

  start_idx <- which(grepl("^\\s*\\d{6}\\s*,", lines))[1]
  if (is.na(start_idx)) stop("Could not find data start in: ", path)

  after <- lines[(start_idx + 1):length(lines)]
  blank_rel <- which(trimws(after) == "")[1]
  end_idx <- if (!is.na(blank_rel)) (start_idx + blank_rel - 1) else length(lines)

  dt <- fread(text = paste(lines[start_idx:end_idx], collapse = "\n"), header = FALSE)
  setDT(dt)

  setnames(dt, 1, "yyyymm")
  dt[, yyyymm := as.integer(yyyymm)]
  dt <- dt[!is.na(yyyymm)]

  dt[, year := yyyymm %/% 100L]
  dt[, mon := yyyymm %% 100L]
  dt[, mdate := as.Date(sprintf("%04d-%02d-01", year, mon))]
  dt[, c("yyyymm", "year", "mon") := NULL]

  ret_cols <- setdiff(names(dt), "mdate")
  for (cc in ret_cols) dt[, (cc) := as.numeric(get(cc)) / 100]

  ports <- c(
    paste0("S", rep(1:4, each = 4), rep(1:4, times = 4)),
    paste0("B", rep(1:4, each = 4), rep(1:4, times = 4))
  )

  if (length(ret_cols) < 32) stop("Expected at least 32 return columns in: ", path)
  setnames(dt, ret_cols[1:32], ports)

  setcolorder(dt, c("mdate", ports))
  dt
}

## Helper: parse generic 5x5 portfolio names ##################################
parse_port_name <- function(x) {
  z <- gsub("[^A-Za-z0-9]", "", x)
  size <- sub("^(S[1-5]).*$", "\\1", z)
  rest <- sub("^S[1-5]", "", z)
  char <- rest
  list(size = size, char = char, port = z)
}

## Helper: wide to long 5x5 ###################################################
wide_to_long_5x5 <- function(wide_dt) {
  stopifnot("mdate" %in% names(wide_dt))

  port_cols <- setdiff(names(wide_dt), "mdate")
  long <- melt(
    wide_dt,
    id.vars = "mdate",
    measure.vars = port_cols,
    variable.name = "port_raw",
    value.name = "ret"
  )

  pp <- lapply(long$port_raw, parse_port_name)
  long[, size := vapply(pp, `[[`, character(1), "size")]
  long[, char := vapply(pp, `[[`, character(1), "char")]
  long[, port := vapply(pp, `[[`, character(1), "port")]

  long
}

## Helper: wide to long 2x4x4 #################################################
wide_to_long_2x4x4_simple <- function(wide_dt) {
  stopifnot("mdate" %in% names(wide_dt))

  port_cols <- setdiff(names(wide_dt), "mdate")
  long <- melt(
    wide_dt,
    id.vars = "mdate",
    measure.vars = port_cols,
    variable.name = "port",
    value.name = "ret"
  )

  long[, size := substr(port, 1, 1)]
  long[, c1 := substr(port, 2, 2)]
  long[, c2 := substr(port, 3, 3)]

  long
}

## Helper: mean matrix from long 5x5 data #####################################
make_mean_matrix <- function(long_dt) {
  tmp <- long_dt[, .(mean_ret = mean(ret, na.rm = TRUE)), by = .(size, char)]
  tmp[, size := factor(size, levels = paste0("S", 1:5))]

  tmp[, char_clean := gsub("[^A-Za-z0-9]", "", as.character(char))]
  tmp[, char_prefix := gsub("[0-9]", "", char_clean)]
  tmp[, char_num := as.integer(gsub("[^0-9]", "", char_clean))]

  pref <- tmp[!is.na(char_prefix) & char_prefix != "", char_prefix][1]
  if (is.na(pref) || pref == "") stop("Could not infer characteristic prefix.")

  tmp[, char := factor(paste0(pref, char_num), levels = paste0(pref, 1:5))]

  mat <- dcast(tmp, size ~ char, value.var = "mean_ret", drop = FALSE)

  needed_cols <- paste0(pref, 1:5)
  for (cc in needed_cols) {
    if (!cc %in% names(mat)) mat[, (cc) := NA_real_]
  }

  mat <- mat[, c("size", needed_cols), with = FALSE]
  out <- as.matrix(mat[, -1])
  rownames(out) <- as.character(mat$size)
  out
}

## Helper: H-L spreads by size ################################################
hl_spreads_by_size <- function(long_dt) {
  tmp <- long_dt[, .(mean_ret = mean(as.numeric(ret), na.rm = TRUE)), by = .(size, char)]
  tmp[, idx := as.integer(gsub("[^0-9]", "", as.character(char)))]
  tmp <- tmp[!is.na(idx)]

  out <- tmp[,
    {
      lo_idx <- min(idx)
      hi_idx <- max(idx)
      lo <- mean_ret[idx == lo_idx][1]
      hi <- mean_ret[idx == hi_idx][1]
      .(low = lo, high = hi, spread = hi - lo)
    },
    by = size
  ]

  out[order(size)]
}

## Helper: time-series correlation by portfolio ###############################
corr_by_port <- function(ff_long, my_long) {
  merged <- merge(
    ff_long[, .(mdate, port, ff_ret = ret)],
    my_long[, .(mdate, port, my_ret = ret)],
    by = c("mdate", "port"),
    all = FALSE
  )

  merged[, .(
    n = .N,
    corr = cor(ff_ret, my_ret, use = "complete.obs"),
    mean_diff = mean(my_ret - ff_ret, na.rm = TRUE)
  ), by = port][order(port)]
}

## Helper: portfolio-level replication regressions ############################
portfolio_replication_regs <- function(ff_long, my_long, lags = NW_LAGS) {
  merged <- merge(
    ff_long[, .(mdate, port, ff_ret = ret)],
    my_long[, .(mdate, port, my_ret = ret)],
    by = c("mdate", "port"),
    all = FALSE
  )

  regs <- lapply(split(merged, by = "port"), function(dt_port) {
    reg <- run_ts_reg_nw(my_ret ~ ff_ret, dt_port[, .(my_ret, ff_ret)], lags = lags)
    if (is.null(reg)) {
      return(NULL)
    }

    data.table(
      port = unique(dt_port$port),
      alpha = get_reg_value(reg, "(Intercept)", "estimate"),
      t_alpha = get_reg_value(reg, "(Intercept)", "tstat"),
      p_alpha = get_reg_value(reg, "(Intercept)", "pval"),
      alpha_stars = get_reg_star(reg, "(Intercept)"),
      beta = get_reg_value(reg, "ff_ret", "estimate"),
      t_beta = get_reg_value(reg, "ff_ret", "tstat"),
      p_beta = get_reg_value(reg, "ff_ret", "pval"),
      beta_stars = get_reg_star(reg, "ff_ret"),
      adj_r2 = get_reg_value(reg, "(Intercept)", "adj_r2"),
      n_obs = nrow(dt_port)
    )
  })

  rbindlist(regs, fill = TRUE)[order(port)]
}

## Helper: 5x5 comparison block ###############################################
compare_block_5x5 <- function(name, ff_long, my_long) {
  cat("\n============================================================\n")
  cat("COMPARISON:", name, "\n")
  cat("Window:", as.character(PORT_START_DATE), "to", as.character(PORT_END_DATE), "\n")
  cat("============================================================\n")

  ff_mat <- make_mean_matrix(ff_long)
  my_mat <- make_mean_matrix(my_long)

  cat("\nFrench mean returns (percent per month)\n")
  print(round(100 * ff_mat, 3))

  cat("\nYour mean returns (percent per month)\n")
  print(round(100 * my_mat, 3))

  cat("\nDifference: Yours minus French (bps per month)\n")
  print(round(10000 * (my_mat - ff_mat), 1))

  cat("\nH-L spreads by size (French, % per month)\n")
  tmpF <- hl_spreads_by_size(ff_long)[, .(size, spread_pct = 100 * spread)]
  print(tmpF[, .(size, spread_pct = round(spread_pct, 3))])

  cat("\nH-L spreads by size (Yours, % per month)\n")
  tmpM <- hl_spreads_by_size(my_long)[, .(size, spread_pct = 100 * spread)]
  print(tmpM[, .(size, spread_pct = round(spread_pct, 3))])

  cat("\nPer-portfolio time-series match\n")
  ts <- corr_by_port(ff_long, my_long)
  print(round_dt(ts, 4))

  cat("\nPer-portfolio replication regressions: Own on French\n")
  reg_tab <- portfolio_replication_regs(ff_long, my_long, lags = NW_LAGS)
  print(round_dt(reg_tab, 4))

  invisible(list(corr = ts, reg = reg_tab))
}

## Helper: 2x4x4 comparison block #############################################
compare_block_2x4x4 <- function(ff_long, my_long, label) {
  merged <- merge(
    ff_long[, .(mdate, port, ff_ret = ret)],
    my_long[, .(mdate, port, my_ret = ret)],
    by = c("mdate", "port"),
    all = FALSE
  )

  cat("\n============================================================\n")
  cat("COMPARISON:", label, "\n")
  cat("Window:", as.character(PORT_START_DATE), "to", as.character(PORT_END_DATE), "\n")
  cat("============================================================\n")

  out <- merged[, .(
    n_obs = .N,
    corr = cor(ff_ret, my_ret, use = "complete.obs"),
    mean_diff = mean(my_ret - ff_ret, na.rm = TRUE)
  ), by = port][order(port)]

  cat("\nPer-portfolio time-series match\n")
  print(round_dt(out, 4))

  cat("\nPer-portfolio replication regressions: Own on French\n")
  reg_tab <- portfolio_replication_regs(ff_long, my_long, lags = NW_LAGS)
  print(round_dt(reg_tab, 4))

  invisible(list(corr = out, reg = reg_tab))
}

## A.1 Comparison to Kenneth French Data: FF5 Factors #########################
kf <- fread(
  file.path(RAW_DATA_DIR, "FF5 Download - 1963_07-2024_12.csv"),
  skip = "mktrf"
)

# The downloaded file is an RTF-wrapped CSV with decimal factor returns and a
# daily date in its final column. Remove the RTF control characters retained in
# the first and last field names/values before standardizing the factor names.
setnames(kf, gsub("[^[:alnum:]_]", "", names(kf)))
setnames(
  kf,
  old = c("f0fs24cf0mktrf", "smb", "hml", "rmw", "cma", "rf", "dateff"),
  new = c("MKT_RF", "SMB", "HML", "RMW", "CMA", "RF", "date_raw"),
  skip_absent = TRUE
)

kf[, date_raw := gsub("[^0-9-]", "", date_raw)]
kf <- kf[grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_raw)]
kf[, mdate := floor_date(as.Date(date_raw), "month")]

cols <- c("MKT_RF", "SMB", "HML", "RMW", "CMA", "RF")
kf[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

kf[, MKT := MKT_RF]
kf <- kf[, .(mdate, MKT, SMB, HML, RMW, CMA)]
kf[, mdate := floor_date(mdate, "month")]

ff5_cmp <- filter_window(copy(ff5), START_DATE, END_DATE)
kf_cmp <- filter_window(copy(kf), START_DATE, END_DATE)

ff_compare <- merge(ff5_cmp, kf_cmp, by = "mdate", all = FALSE)

setnames(
  ff_compare,
  old = c(
    "SMB.x", "HML.x", "RMW.x", "CMA.x",
    "SMB.y", "HML.y", "RMW.y", "CMA.y"
  ),
  new = c(
    "SMB.own", "HML.own", "RMW.own", "CMA.own",
    "SMB.kf", "HML.kf", "RMW.kf", "CMA.kf"
  )
)

## A.1.1 Mean Returns #########################################################
mean_comp <-
  ff_compare[, .(
    MKT_own = mean(MKT.x, na.rm = TRUE) * 100,
    MKT_kf  = mean(MKT.y, na.rm = TRUE) * 100,
    SMB_own = mean(SMB.own, na.rm = TRUE) * 100,
    SMB_kf  = mean(SMB.kf, na.rm = TRUE) * 100,
    HML_own = mean(HML.own, na.rm = TRUE) * 100,
    HML_kf  = mean(HML.kf, na.rm = TRUE) * 100,
    RMW_own = mean(RMW.own, na.rm = TRUE) * 100,
    RMW_kf  = mean(RMW.kf, na.rm = TRUE) * 100,
    CMA_own = mean(CMA.own, na.rm = TRUE) * 100,
    CMA_kf  = mean(CMA.kf, na.rm = TRUE) * 100
  )]

mean_long <-
  melt(mean_comp, variable.name = "series", value.name = "mean_pct")

mean_long[, source := fifelse(grepl("_own$", series), "Own", "KF")]
mean_long[, factor := gsub("_(own|kf)$", "", series)]
mean_long[, series := NULL]

mean_diff <-
  dcast(mean_long, factor ~ source, value.var = "mean_pct")

mean_diff[, Difference := Own - KF]

cat("\n============================================================\n")
cat("A.1.1 FF5 Mean Return Comparison (% per month)\n")
cat("============================================================\n")
print(mean_diff)

## A.1.2 Correlations #########################################################
cor_mat <-
  round(
    cor(
      ff_compare[, .(
        MKT.x,
        SMB.own, HML.own, RMW.own, CMA.own,
        MKT.y,
        SMB.kf,  HML.kf,  RMW.kf,  CMA.kf
      )],
      use = "pairwise.complete.obs"
    ),
    3
  )

ff5_cor_compare <-
  as.data.table(
    cor_mat[
      c("MKT.x", "SMB.own", "HML.own", "RMW.own", "CMA.own"),
      c("MKT.y", "SMB.kf", "HML.kf", "RMW.kf", "CMA.kf")
    ],
    keep.rownames = "own_factor"
  )

cat("\n============================================================\n")
cat("A.1.2 FF5 Correlations: Own vs Kenneth French\n")
cat("============================================================\n")
print(ff5_cor_compare)

## A.1.3 Time-Series Replication Regressions ##################################
reg_MKT <- run_ts_reg_nw(MKT.x ~ MKT.y, ff_compare[, .(MKT.x, MKT.y)], NW_LAGS)
reg_SMB <- run_ts_reg_nw(SMB.own ~ SMB.kf, ff_compare[, .(SMB.own, SMB.kf)], NW_LAGS)
reg_HML <- run_ts_reg_nw(HML.own ~ HML.kf, ff_compare[, .(HML.own, HML.kf)], NW_LAGS)
reg_RMW <- run_ts_reg_nw(RMW.own ~ RMW.kf, ff_compare[, .(RMW.own, RMW.kf)], NW_LAGS)
reg_CMA <- run_ts_reg_nw(CMA.own ~ CMA.kf, ff_compare[, .(CMA.own, CMA.kf)], NW_LAGS)

ff5_regression_table <- rbindlist(list(
  cbind(data.table(factor = "MKT"), reg_MKT),
  cbind(data.table(factor = "SMB"), reg_SMB),
  cbind(data.table(factor = "HML"), reg_HML),
  cbind(data.table(factor = "RMW"), reg_RMW),
  cbind(data.table(factor = "CMA"), reg_CMA)
), fill = TRUE)

cat("\n============================================================\n")
cat("A.1.3 FF5 Replication Regressions (Own on Kenneth French)\n")
cat("============================================================\n")
print(round_dt(ff5_regression_table, 4))

## A.1.4 Compact Replication Summary ##########################################
ff5_regression_compact <-
  data.table(
    factor = c("MKT", "SMB", "HML", "RMW", "CMA"),
    alpha = c(
      get_reg_value(reg_MKT, "(Intercept)", "estimate"),
      get_reg_value(reg_SMB, "(Intercept)", "estimate"),
      get_reg_value(reg_HML, "(Intercept)", "estimate"),
      get_reg_value(reg_RMW, "(Intercept)", "estimate"),
      get_reg_value(reg_CMA, "(Intercept)", "estimate")
    ),
    t_alpha = c(
      get_reg_value(reg_MKT, "(Intercept)", "tstat"),
      get_reg_value(reg_SMB, "(Intercept)", "tstat"),
      get_reg_value(reg_HML, "(Intercept)", "tstat"),
      get_reg_value(reg_RMW, "(Intercept)", "tstat"),
      get_reg_value(reg_CMA, "(Intercept)", "tstat")
    ),
    p_alpha = c(
      get_reg_value(reg_MKT, "(Intercept)", "pval"),
      get_reg_value(reg_SMB, "(Intercept)", "pval"),
      get_reg_value(reg_HML, "(Intercept)", "pval"),
      get_reg_value(reg_RMW, "(Intercept)", "pval"),
      get_reg_value(reg_CMA, "(Intercept)", "pval")
    ),
    alpha_stars = c(
      get_reg_star(reg_MKT, "(Intercept)"),
      get_reg_star(reg_SMB, "(Intercept)"),
      get_reg_star(reg_HML, "(Intercept)"),
      get_reg_star(reg_RMW, "(Intercept)"),
      get_reg_star(reg_CMA, "(Intercept)")
    ),
    beta = c(
      get_reg_value(reg_MKT, "MKT.y", "estimate"),
      get_reg_value(reg_SMB, "SMB.kf", "estimate"),
      get_reg_value(reg_HML, "HML.kf", "estimate"),
      get_reg_value(reg_RMW, "RMW.kf", "estimate"),
      get_reg_value(reg_CMA, "CMA.kf", "estimate")
    ),
    t_beta = c(
      get_reg_value(reg_MKT, "MKT.y", "tstat"),
      get_reg_value(reg_SMB, "SMB.kf", "tstat"),
      get_reg_value(reg_HML, "HML.kf", "tstat"),
      get_reg_value(reg_RMW, "RMW.kf", "tstat"),
      get_reg_value(reg_CMA, "CMA.kf", "tstat")
    ),
    p_beta = c(
      get_reg_value(reg_MKT, "MKT.y", "pval"),
      get_reg_value(reg_SMB, "SMB.kf", "pval"),
      get_reg_value(reg_HML, "HML.kf", "pval"),
      get_reg_value(reg_RMW, "RMW.kf", "pval"),
      get_reg_value(reg_CMA, "CMA.kf", "pval")
    ),
    beta_stars = c(
      get_reg_star(reg_MKT, "MKT.y"),
      get_reg_star(reg_SMB, "SMB.kf"),
      get_reg_star(reg_HML, "HML.kf"),
      get_reg_star(reg_RMW, "RMW.kf"),
      get_reg_star(reg_CMA, "CMA.kf")
    ),
    adj_r2 = c(
      get_reg_value(reg_MKT, "(Intercept)", "adj_r2"),
      get_reg_value(reg_SMB, "(Intercept)", "adj_r2"),
      get_reg_value(reg_HML, "(Intercept)", "adj_r2"),
      get_reg_value(reg_RMW, "(Intercept)", "adj_r2"),
      get_reg_value(reg_CMA, "(Intercept)", "adj_r2")
    )
  )

cat("\n============================================================\n")
cat("A.1.4 FF5 Compact Replication Summary\n")
cat("============================================================\n")
print(round_dt(ff5_regression_compact, 4))

## A.2 5x5 Comparison to Kenneth French Data ##################################
ff_bm_file <- file.path(RAW_DATA_DIR, "kf_portfolios", "25_Portfolios_5x5.csv")
ff_op_file <- file.path(RAW_DATA_DIR, "kf_portfolios", "25_Portfolios_ME_OP_5x5.csv")
ff_inv_file <- file.path(RAW_DATA_DIR, "kf_portfolios", "25_Portfolios_ME_INV_5x5.csv")

my_bm_rds <- file.path(APPENDIX_DIR, "my_Size_BM_5x5_wide.rds")
my_op_rds <- file.path(APPENDIX_DIR, "my_Size_OP_5x5_wide.rds")
my_inv_rds <- file.path(APPENDIX_DIR, "my_Size_INV_5x5_wide.rds")

ff_bm <- read_french_5x5(ff_bm_file, char_prefix = "V")
ff_op <- read_french_5x5(ff_op_file, char_prefix = "P")
ff_inv <- read_french_5x5(ff_inv_file, char_prefix = "I")

my_bm_w <- readRDS(my_bm_rds)
my_op_w <- readRDS(my_op_rds)
my_inv_w <- readRDS(my_inv_rds)

setDT(my_bm_w)
setDT(my_op_w)
setDT(my_inv_w)

my_bm_w[, mdate := as.Date(mdate)]
my_op_w[, mdate := as.Date(mdate)]
my_inv_w[, mdate := as.Date(mdate)]

ff_bm_L <- wide_to_long_5x5(ff_bm)
ff_op_L <- wide_to_long_5x5(ff_op)
ff_inv_L <- wide_to_long_5x5(ff_inv)

my_bm_L <- wide_to_long_5x5(my_bm_w)
my_op_L <- wide_to_long_5x5(my_op_w)
my_inv_L <- wide_to_long_5x5(my_inv_w)

ff_bm_L <- filter_window(ff_bm_L, PORT_START_DATE, PORT_END_DATE)
ff_op_L <- filter_window(ff_op_L, PORT_START_DATE, PORT_END_DATE)
ff_inv_L <- filter_window(ff_inv_L, PORT_START_DATE, PORT_END_DATE)

my_bm_L <- filter_window(my_bm_L, PORT_START_DATE, PORT_END_DATE)
my_op_L <- filter_window(my_op_L, PORT_START_DATE, PORT_END_DATE)
my_inv_L <- filter_window(my_inv_L, PORT_START_DATE, PORT_END_DATE)

bm_5x5_results <- compare_block_5x5("Size x B/M 5x5", ff_bm_L, my_bm_L)
op_5x5_results <- compare_block_5x5("Size x OP 5x5", ff_op_L, my_op_L)
inv_5x5_results <- compare_block_5x5("Size x INV 5x5", ff_inv_L, my_inv_L)

## A.3 2x4x4 Comparison to Kenneth French Data ################################
ff_op_inv <- read_french_2x4x4(file.path(RAW_DATA_DIR, "kf_portfolios", "32_Portfolios_ME_OP_INV_2x4x4.csv"))
ff_bm_inv <- read_french_2x4x4(file.path(RAW_DATA_DIR, "kf_portfolios", "32_Portfolios_ME_BEME_INV_2x4x4.csv"))
ff_bm_op <- read_french_2x4x4(file.path(RAW_DATA_DIR, "kf_portfolios", "32_Portfolios_ME_BEME_OP_2x4x4.csv"))

ff_op_inv_L <- wide_to_long_2x4x4_simple(ff_op_inv)
ff_bm_inv_L <- wide_to_long_2x4x4_simple(ff_bm_inv)
ff_bm_op_L <- wide_to_long_2x4x4_simple(ff_bm_op)

my_op_inv_w <- readRDS(file.path(APPENDIX_DIR, "my_Size_OP_INV_2x4x4_wide.rds"))
my_bm_inv_w <- readRDS(file.path(APPENDIX_DIR, "my_Size_BM_INV_2x4x4_wide.rds"))
my_bm_op_w <- readRDS(file.path(APPENDIX_DIR, "my_Size_BM_OP_2x4x4_wide.rds"))

setDT(my_op_inv_w)
setDT(my_bm_inv_w)
setDT(my_bm_op_w)

my_op_inv_w[, mdate := as.Date(mdate)]
my_bm_inv_w[, mdate := as.Date(mdate)]
my_bm_op_w[, mdate := as.Date(mdate)]

my_op_inv_L <- wide_to_long_2x4x4_simple(my_op_inv_w)
my_bm_inv_L <- wide_to_long_2x4x4_simple(my_bm_inv_w)
my_bm_op_L <- wide_to_long_2x4x4_simple(my_bm_op_w)

ff_op_inv_L <- filter_window(ff_op_inv_L, PORT_START_DATE, PORT_END_DATE)
ff_bm_inv_L <- filter_window(ff_bm_inv_L, PORT_START_DATE, PORT_END_DATE)
ff_bm_op_L <- filter_window(ff_bm_op_L, PORT_START_DATE, PORT_END_DATE)

my_op_inv_L <- filter_window(my_op_inv_L, PORT_START_DATE, PORT_END_DATE)
my_bm_inv_L <- filter_window(my_bm_inv_L, PORT_START_DATE, PORT_END_DATE)
my_bm_op_L <- filter_window(my_bm_op_L, PORT_START_DATE, PORT_END_DATE)

op_inv_2x4x4_results <- compare_block_2x4x4(ff_op_inv_L, my_op_inv_L, "Size x OP x INV")
bm_inv_2x4x4_results <- compare_block_2x4x4(ff_bm_inv_L, my_bm_inv_L, "Size x BM x INV")
bm_op_2x4x4_results <- compare_block_2x4x4(ff_bm_op_L, my_bm_op_L, "Size x BM x OP")

## A.4 Factor Diagnostics #####################################################
cat("\n============================================================\n")
cat("A.4.1 Factor Date Ranges\n")
cat("============================================================\n")
print(range(ff5$mdate))
print(range(SMQ$mdate))

cat("\n============================================================\n")
cat("A.4.2 FF5 Summary Statistics\n")
cat("============================================================\n")
print(summary(ff5[, .(MKT, SMB, HML, RMW, CMA)]))
print(apply(ff5[, .(MKT, SMB, HML, RMW, CMA)], 2, sd, na.rm = TRUE))

cat("\n============================================================\n")
cat("A.4.3 SMQ Summary Statistics\n")
cat("============================================================\n")
print(summary(SMQ[, -1]))
print(apply(SMQ[, -1], 2, sd, na.rm = TRUE))

par(mfrow = c(3, 2), mar = c(3, 4, 2, 1))
plot(ff5$mdate, ff5$MKT, type = "l", main = "MKT")
plot(ff5$mdate, ff5$SMB, type = "l", main = "SMB")
plot(ff5$mdate, ff5$HML, type = "l", main = "HML")
plot(ff5$mdate, ff5$RMW, type = "l", main = "RMW")
plot(ff5$mdate, ff5$CMA, type = "l", main = "CMA")
par(mfrow = c(1, 1))

par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(SMQ$mdate, SMQ$SMQ_firm_VW,
  type = "l", lwd = 1.5,
  main = "SMQ - Firm Level (VW)", ylab = "Return", xlab = ""
)
plot(SMQ$mdate, SMQ$SMQ_sic_VW,
  type = "l", lwd = 1.5,
  main = "SMQ - SIC Level (VW)", ylab = "Return", xlab = ""
)
plot(SMQ$mdate, SMQ$SMQ_sic2_VW,
  type = "l", lwd = 1.5,
  main = "SMQ - SIC2 Level (VW)", ylab = "Return", xlab = "Time"
)
par(mfrow = c(1, 1))

par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(SMQ$mdate, SMQ$SMQ_firm_K,
  type = "l", lwd = 1.5,
  main = "SMQ - Firm Level (K-weighted)", ylab = "Return", xlab = ""
)
plot(SMQ$mdate, SMQ$SMQ_sic_K,
  type = "l", lwd = 1.5,
  main = "SMQ - SIC Level (K-weighted)", ylab = "Return", xlab = ""
)
plot(SMQ$mdate, SMQ$SMQ_sic2_K,
  type = "l", lwd = 1.5,
  main = "SMQ - SIC2 Level (K-weighted)", ylab = "Return", xlab = "Time"
)
par(mfrow = c(1, 1))

par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(SMQ$mdate, cum_safe(SMQ$SMQ_firm_VW),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - Firm Level (VW)", ylab = "Cumulative Return", xlab = ""
)
plot(SMQ$mdate, cum_safe(SMQ$SMQ_sic_VW),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - SIC Level (VW)", ylab = "Cumulative Return", xlab = ""
)
plot(SMQ$mdate, cum_safe(SMQ$SMQ_sic2_VW),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - SIC2 Level (VW)", ylab = "Cumulative Return", xlab = "Time"
)
par(mfrow = c(1, 1))

par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(SMQ$mdate, cum_safe(SMQ$SMQ_firm_K),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - Firm Level (K-weighted)", ylab = "Cumulative Return", xlab = ""
)
plot(SMQ$mdate, cum_safe(SMQ$SMQ_sic_K),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - SIC Level (K-weighted)", ylab = "Cumulative Return", xlab = ""
)
plot(SMQ$mdate, cum_safe(SMQ$SMQ_sic2_K),
  type = "l", lwd = 1.5,
  main = "Cumulative SMQ - SIC2 Level (K-weighted)", ylab = "Cumulative Return", xlab = "Time"
)
par(mfrow = c(1, 1))

## A.5 Concentration Diagnostics ##############################################
## A.5.1 Number of Firms per SIC ##############################################
cat("\n============================================================\n")
cat("A.5.1 Number of Firms per SIC\n")
cat("============================================================\n")
print(
  data[!is.na(sic), .N, by = sic][
    , .(
      mean = mean(N),
      p10 = quantile(N, 0.10),
      median = median(N),
      p90 = quantile(N, 0.90),
      max = max(N)
    )
  ]
)

## A.5.2 Number of Firms per SIC2 #############################################
cat("\n============================================================\n")
cat("A.5.2 Number of Firms per SIC2\n")
cat("============================================================\n")
print(
  data[!is.na(sic2), .N, by = sic2][
    , .(
      mean = mean(N),
      p10 = quantile(N, 0.10),
      median = median(N),
      p90 = quantile(N, 0.90),
      max = max(N)
    )
  ]
)

## A.5.3 Top SIC Shares by Market Equity ######################################
cat("\n============================================================\n")
cat("A.5.3 Top SIC Shares by Market Equity\n")
cat("============================================================\n")
print(
  data[
    month == 7 & !is.na(me_clean) & !is.na(sic),
    .(ME = sum(me_clean, na.rm = TRUE)),
    by = sic
  ][
    , share := ME / sum(ME)
  ][order(-share)][1:10]
)

## A.5.4 Top SIC Shares by Capital ############################################
cat("\n============================================================\n")
cat("A.5.4 Top SIC Shares by Capital\n")
cat("============================================================\n")
print(
  data[
    month == 7 & !is.na(ppegt) & !is.na(sic),
    .(K = sum(ppegt, na.rm = TRUE)),
    by = sic
  ][
    , share := K / sum(K)
  ][order(-share)][1:10]
)

## A.5.5 Average Firms per SIC over Time ######################################
cat("\n============================================================\n")
cat("A.5.5 Average Firms per SIC over Time\n")
cat("============================================================\n")
print(
  data[!is.na(sic), .N, by = .(year, sic)][
    , .(mean_firms = mean(N)),
    by = year
  ]
)

## A.5.6 Average Firms per SIC2 over Time #####################################
cat("\n============================================================\n")
cat("A.5.6 Average Firms per SIC2 over Time\n")
cat("============================================================\n")
print(
  data[!is.na(sic2), .N, by = .(year, sic2)][
    , .(mean_firms = mean(N)),
    by = year
  ]
)

## A.6 Characteristic Boxplot Diagnostics #####################################
form_sample <- data[
  month == 6 &
    has_dec == TRUE &
    has_june == TRUE &
    mdate >= START_DATE &
    mdate <= END_DATE
]

form_sample[, decade :=
  fifelse(
    ffyear < 1990, "1980s",
    fifelse(
      ffyear < 2000, "1990s",
      fifelse(
        ffyear < 2010, "2000s",
        fifelse(ffyear < 2020, "2010s", "2020s")
      )
    )
  )]

form_sample[, decade := factor(
  decade,
  levels = c("1980s", "1990s", "2000s", "2010s", "2020s")
)]

plot_vars <- c("bm_ff", "op_clean", "inv_annual", "chi_firm", "chi_sic", "chi_sic2")

for (v in plot_vars) {
  form_sample[, paste0(v, "_w1") := winsor_1pct(get(v))]
}

par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))
for (v in plot_vars) {
  boxplot(
    form_sample[[v]] ~ form_sample$decade,
    main = paste(v, "by Decade"),
    xlab = "",
    ylab = v
  )
}
par(mfrow = c(1, 1))

par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))
for (v in plot_vars) {
  vw <- paste0(v, "_w1")

  boxplot(
    form_sample[[vw]] ~ form_sample$decade,
    main = paste(v, "by Decade (winsorized 1%)"),
    xlab = "",
    ylab = v
  )
}
par(mfrow = c(1, 1))

## Optional Saves #############################################################
saveRDS(ff5_cor_compare, "appendix_ff5_cor_compare.rds")
saveRDS(ff5_regression_table, "appendix_ff5_regression_table.rds")
saveRDS(ff5_regression_compact, "appendix_ff5_regression_compact.rds")

saveRDS(bm_5x5_results, "appendix_bm_5x5_results.rds")
saveRDS(op_5x5_results, "appendix_op_5x5_results.rds")
saveRDS(inv_5x5_results, "appendix_inv_5x5_results.rds")

saveRDS(op_inv_2x4x4_results, "appendix_op_inv_2x4x4_results.rds")
saveRDS(bm_inv_2x4x4_results, "appendix_bm_inv_2x4x4_results.rds")
saveRDS(bm_op_2x4x4_results, "appendix_bm_op_2x4x4_results.rds")







## A.7 Examine Negative Mean Excess Return Portfolio ###########################
## Further examination of S14 = Small firms, OP quartile 1, Investment quartile 4
## This porftolio earns a mean excess return of -0.002
### Examine S14 Return Dynamics and Adjustment Costs ###########################

library(data.table)
library(lubridate)
library(ggplot2)

### Settings ##################################################################

START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")

PORT_COL <- "port_OP_INV_2x4x4"
NEG_PORT_LABEL <- "S14"

RETURNS_ARE_EXCESS <- FALSE

### Helper Functions ###########################################################

mean_safe <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

median_safe <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

sd_safe <- function(x) {
  if (sum(!is.na(x)) < 2) {
    return(NA_real_)
  }
  sd(x, na.rm = TRUE)
}

period_return_safe <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  prod(1 + x) - 1
}

wmean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

winsor_vec <- function(x, p = 0.01) {
  if (sum(!is.na(x)) < 20) {
    return(x)
  }
  q <- quantile(x, probs = c(p, 1 - p), na.rm = TRUE, type = 7)
  pmin(pmax(x, q[1]), q[2])
}

### Load Firm-Level Data and Risk-Free Rate ###################################

panel_labels <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(panel_labels)

panel_labels[, mdate := floor_date(as.Date(mdate), "month")]

rf_monthly <- panel_labels[
  mdate >= START_DATE &
    mdate <= END_DATE &
    !is.na(RF),
  .(RF = mean(RF, na.rm = TRUE)),
  by = mdate
]

### Load 32 Portfolio Returns #################################################

port_wide <- readRDS("my_Size_OP_INV_2x4x4_wide.rds")
setDT(port_wide)

port_wide[, mdate := floor_date(as.Date(mdate), "month")]

port_cols <- grep("^[SB][1-4][1-4]$", names(port_wide), value = TRUE)

if (!NEG_PORT_LABEL %in% port_cols) {
  stop("S14 is not found in my_Size_OP_INV_2x4x4_wide.rds.")
}

ret_long <- melt(
  port_wide[
    mdate >= START_DATE &
      mdate <= END_DATE,
    c("mdate", port_cols),
    with = FALSE
  ],
  id.vars = "mdate",
  variable.name = "port",
  value.name = "port_ret"
)

ret_long[, port := as.character(port)]
ret_long[, size := substr(port, 1, 1)]
ret_long[, op_q := as.integer(substr(port, 2, 2))]
ret_long[, inv_q := as.integer(substr(port, 3, 3))]

ret_long <- merge(
  ret_long,
  rf_monthly,
  by = "mdate",
  all.x = TRUE
)

ret_long[, excess_ret_minus_rf := port_ret - RF]

if (RETURNS_ARE_EXCESS) {
  ret_long[, excess_ret := port_ret]
} else {
  ret_long[, excess_ret := port_ret - RF]
}

### Check Whether the File Contains Raw or Excess Returns ######################

s14_return_check <- ret_long[
  port == NEG_PORT_LABEL,
  .(
    mean_port_ret_pct = 100 * mean_safe(port_ret),
    mean_port_ret_minus_rf_pct = 100 * mean_safe(excess_ret_minus_rf),
    mean_rf_pct = 100 * mean_safe(RF),
    n_months = .N
  )
]

cat("\nS14 return check:\n")
print(s14_return_check)

cat("\nIf mean_port_ret_pct matches Table 7, set RETURNS_ARE_EXCESS = TRUE.\n")
cat("If mean_port_ret_minus_rf_pct matches Table 7, keep RETURNS_ARE_EXCESS = FALSE.\n")

### Return Summary Across All 32 Portfolios ###################################

ret_summary <- ret_long[
  ,
  .(
    mean_excess_pct = 100 * mean_safe(excess_ret),
    median_excess_pct = 100 * median_safe(excess_ret),
    sd_excess_pct = 100 * sd_safe(excess_ret),
    min_excess_pct = 100 * min(excess_ret, na.rm = TRUE),
    max_excess_pct = 100 * max(excess_ret, na.rm = TRUE),
    cumulative_excess_pct = 100 * period_return_safe(excess_ret),
    n_months = .N
  ),
  by = .(port, size, op_q, inv_q)
][order(mean_excess_pct)]

cat("\nAll 32 portfolios ranked by mean excess return:\n")
print(ret_summary)

cat("\nS14 relative return position:\n")
print(ret_summary[port == NEG_PORT_LABEL])

### S14 Return Dynamics Over Time #############################################

s14_ts <- ret_long[port == NEG_PORT_LABEL][order(mdate)]

s14_ts[, excess_ret_pct := 100 * excess_ret]
s14_ts[, cum_excess_ret_pct := 100 * (cumprod(1 + fifelse(is.na(excess_ret), 0, excess_ret)) - 1)]

s14_ts[, roll_12_mean_pct := 100 * frollmean(excess_ret, n = 12, align = "right", na.rm = TRUE)]
s14_ts[, roll_36_mean_pct := 100 * frollmean(excess_ret, n = 36, align = "right", na.rm = TRUE)]
s14_ts[, roll_60_mean_pct := 100 * frollmean(excess_ret, n = 60, align = "right", na.rm = TRUE)]

s14_ts[, year := year(mdate)]

s14_ts[, period :=
  fifelse(
    year <= 1989, "1982-1989",
    fifelse(
      year <= 1999, "1990-1999",
      fifelse(
        year <= 2009, "2000-2009",
        fifelse(year <= 2019, "2010-2019", "2020-2024")
      )
    )
  )]

s14_by_period <- s14_ts[
  ,
  .(
    mean_excess_pct = 100 * mean_safe(excess_ret),
    median_excess_pct = 100 * median_safe(excess_ret),
    sd_excess_pct = 100 * sd_safe(excess_ret),
    cumulative_excess_pct = 100 * period_return_safe(excess_ret),
    min_excess_pct = 100 * min(excess_ret, na.rm = TRUE),
    max_excess_pct = 100 * max(excess_ret, na.rm = TRUE),
    n_months = .N
  ),
  by = period
]

s14_by_year <- s14_ts[
  ,
  .(
    mean_excess_pct = 100 * mean_safe(excess_ret),
    cumulative_excess_pct = 100 * period_return_safe(excess_ret),
    sd_excess_pct = 100 * sd_safe(excess_ret),
    n_months = .N
  ),
  by = year
][order(year)]

cat("\nS14 return by period:\n")
print(s14_by_period)

cat("\nS14 return by year:\n")
print(s14_by_year)

### Plot S14 Returns Over Time ################################################

p_s14_monthly <- ggplot(s14_ts, aes(x = mdate, y = excess_ret_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.4) +
  labs(
    title = "S14 Monthly Excess Returns",
    subtitle = "Small firms, low operating profitability, high investment",
    x = "",
    y = "Monthly excess return (%)"
  ) +
  theme_minimal()

p_s14_rolling <- ggplot(s14_ts, aes(x = mdate)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(aes(y = roll_12_mean_pct), linewidth = 0.4) +
  geom_line(aes(y = roll_36_mean_pct), linewidth = 0.6) +
  geom_line(aes(y = roll_60_mean_pct), linewidth = 0.8) +
  labs(
    title = "S14 Rolling Mean Excess Returns",
    subtitle = "12-, 36-, and 60-month rolling averages",
    x = "",
    y = "Rolling mean excess return (%)"
  ) +
  theme_minimal()

p_s14_cumulative <- ggplot(s14_ts, aes(x = mdate, y = cum_excess_ret_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  labs(
    title = "S14 Cumulative Excess Return",
    subtitle = "Small firms, low operating profitability, high investment",
    x = "",
    y = "Cumulative excess return (%)"
  ) +
  theme_minimal()

print(p_s14_monthly)
print(p_s14_rolling)
print(p_s14_cumulative)

### Adjustment Cost Comparison Across 32 Portfolios ############################

ME_COL <- "me_clean"
if (!ME_COL %in% names(panel_labels)) ME_COL <- "me_ff"

panel_chi <- panel_labels[
  mdate >= START_DATE &
    mdate <= END_DATE &
    !is.na(get(PORT_COL))
]

panel_chi[, port := as.character(get(PORT_COL))]
panel_chi[, size := substr(port, 1, 1)]
panel_chi[, op_q := as.integer(substr(port, 2, 2))]
panel_chi[, inv_q := as.integer(substr(port, 3, 3))]

panel_chi[, chi_firm_w := winsor_vec(chi_firm), by = mdate]
panel_chi[, chi_sic_w := winsor_vec(chi_sic), by = mdate]
panel_chi[, chi_sic2_w := winsor_vec(chi_sic2), by = mdate]

chi_monthly <- panel_chi[
  ,
  .(
    n_firms = uniqueN(gvkey),
    share_chi_firm_available = mean(!is.na(chi_firm)),
    share_chi_sic_available = mean(!is.na(chi_sic)),
    share_chi_sic2_available = mean(!is.na(chi_sic2)),
    chi_firm_ew = mean_safe(chi_firm_w),
    chi_firm_vw = wmean_safe(chi_firm_w, get(ME_COL)),
    chi_firm_median = median_safe(chi_firm),
    chi_sic_ew = mean_safe(chi_sic_w),
    chi_sic_vw = wmean_safe(chi_sic_w, get(ME_COL)),
    chi_sic_median = median_safe(chi_sic),
    chi_sic2_ew = mean_safe(chi_sic2_w),
    chi_sic2_vw = wmean_safe(chi_sic2_w, get(ME_COL)),
    chi_sic2_median = median_safe(chi_sic2),
    mean_op = mean_safe(op_clean),
    mean_inv = mean_safe(inv_annual),
    mean_bm = mean_safe(bm_ff),
    mean_me = mean_safe(get(ME_COL))
  ),
  by = .(mdate, port, size, op_q, inv_q)
]

chi_summary <- chi_monthly[
  ,
  .(
    n_months = .N,
    avg_n_firms = mean_safe(n_firms),
    avg_share_chi_firm_available = mean_safe(share_chi_firm_available),
    avg_share_chi_sic_available = mean_safe(share_chi_sic_available),
    avg_share_chi_sic2_available = mean_safe(share_chi_sic2_available),
    avg_chi_firm_ew = mean_safe(chi_firm_ew),
    avg_chi_firm_vw = mean_safe(chi_firm_vw),
    avg_chi_firm_median = mean_safe(chi_firm_median),
    avg_chi_sic_ew = mean_safe(chi_sic_ew),
    avg_chi_sic_vw = mean_safe(chi_sic_vw),
    avg_chi_sic_median = mean_safe(chi_sic_median),
    avg_chi_sic2_ew = mean_safe(chi_sic2_ew),
    avg_chi_sic2_vw = mean_safe(chi_sic2_vw),
    avg_chi_sic2_median = mean_safe(chi_sic2_median),
    avg_op = mean_safe(mean_op),
    avg_inv = mean_safe(mean_inv),
    avg_bm = mean_safe(mean_bm),
    avg_me = mean_safe(mean_me)
  ),
  by = .(port, size, op_q, inv_q)
]

chi_summary[, rank_high_chi_sic_vw := frank(-avg_chi_sic_vw, ties.method = "average", na.last = "keep")]
chi_summary[, rank_high_chi_sic_ew := frank(-avg_chi_sic_ew, ties.method = "average", na.last = "keep")]
chi_summary[, rank_high_chi_sic2_vw := frank(-avg_chi_sic2_vw, ties.method = "average", na.last = "keep")]
chi_summary[, rank_high_chi_firm_vw := frank(-avg_chi_firm_vw, ties.method = "average", na.last = "keep")]

chi_return_summary <- merge(
  chi_summary,
  ret_summary[, .(port, mean_excess_pct, median_excess_pct, sd_excess_pct, cumulative_excess_pct)],
  by = "port",
  all.x = TRUE
)

cat("\nS14 adjustment cost and return summary:\n")
print(chi_return_summary[port == NEG_PORT_LABEL])

cat("\nAll 32 portfolios ranked by value-weighted SIC-level chi:\n")
print(
  chi_return_summary[
    order(-avg_chi_sic_vw),
    .(
      port,
      size,
      op_q,
      inv_q,
      mean_excess_pct,
      avg_chi_sic_vw,
      rank_high_chi_sic_vw,
      avg_chi_sic_ew,
      avg_chi_sic2_vw,
      avg_chi_firm_vw,
      avg_share_chi_firm_available,
      avg_n_firms
    )
  ]
)

cat("\nSmall portfolios only, ranked by value-weighted SIC-level chi:\n")
print(
  chi_return_summary[
    size == "S"
  ][
    order(-avg_chi_sic_vw),
    .(
      port,
      op_q,
      inv_q,
      mean_excess_pct,
      avg_chi_sic_vw,
      rank_high_chi_sic_vw,
      avg_chi_sic2_vw,
      avg_chi_firm_vw,
      avg_n_firms
    )
  ]
)

cat("\nLow OP portfolios only:\n")
print(
  chi_return_summary[
    op_q == 1
  ][
    order(size, inv_q),
    .(
      port,
      size,
      inv_q,
      mean_excess_pct,
      avg_chi_sic_vw,
      avg_chi_sic2_vw,
      avg_chi_firm_vw,
      avg_n_firms
    )
  ]
)

cat("\nHigh investment portfolios only:\n")
print(
  chi_return_summary[
    inv_q == 4
  ][
    order(size, op_q),
    .(
      port,
      size,
      op_q,
      mean_excess_pct,
      avg_chi_sic_vw,
      avg_chi_sic2_vw,
      avg_chi_firm_vw,
      avg_n_firms
    )
  ]
)

### Save Outputs ##############################################################

saveRDS(s14_ts, "S14_return_timeseries.rds")
saveRDS(s14_by_period, "S14_return_by_period.rds")
saveRDS(s14_by_year, "S14_return_by_year.rds")
saveRDS(ret_summary, "OP_INV_2x4x4_return_summary.rds")
saveRDS(chi_summary, "OP_INV_2x4x4_chi_summary.rds")
saveRDS(chi_return_summary, "OP_INV_2x4x4_chi_return_summary.rds")

### A.7.1 Book-to-Market Rank and Industry Composition #########################

library(data.table)
library(lubridate)

### Settings ##################################################################

START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")

PORT_COL <- "port_OP_INV_2x4x4"
NEG_PORT_LABEL <- "S14"

ME_COL <- "me_clean"

### Helper Functions ###########################################################

mean_safe <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

median_safe <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  median(x, na.rm = TRUE)
}

wmean_safe <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

### Load Data #################################################################

panel_labels <- readRDS("crsp_compustat_with_portfolio_labels.rds")
setDT(panel_labels)

panel_labels[, mdate := floor_date(as.Date(mdate), "month")]

if (!ME_COL %in% names(panel_labels)) ME_COL <- "me_ff"

panel_32 <- panel_labels[
  mdate >= START_DATE &
    mdate <= END_DATE &
    !is.na(get(PORT_COL))
]

panel_32[, port := as.character(get(PORT_COL))]
panel_32[, size := substr(port, 1, 1)]
panel_32[, op_q := as.integer(substr(port, 2, 2))]
panel_32[, inv_q := as.integer(substr(port, 3, 3))]

### Book-to-Market Summary Across 32 Portfolios ###############################

bm_monthly <- panel_32[
  ,
  .(
    n_firms = uniqueN(gvkey),
    bm_ew = mean_safe(bm_ff),
    bm_vw = wmean_safe(bm_ff, get(ME_COL)),
    bm_median = median_safe(bm_ff),
    be_ew = mean_safe(be_clean),
    be_vw = wmean_safe(be_clean, get(ME_COL)),
    be_median = median_safe(be_clean)
  ),
  by = .(mdate, port, size, op_q, inv_q)
]

bm_summary <- bm_monthly[
  ,
  .(
    n_months = .N,
    avg_n_firms = mean_safe(n_firms),
    avg_bm_ew = mean_safe(bm_ew),
    avg_bm_vw = mean_safe(bm_vw),
    avg_bm_median = mean_safe(bm_median),
    avg_be_ew = mean_safe(be_ew),
    avg_be_vw = mean_safe(be_vw),
    avg_be_median = mean_safe(be_median)
  ),
  by = .(port, size, op_q, inv_q)
]

bm_summary[, rank_high_bm_ew := frank(-avg_bm_ew, ties.method = "average", na.last = "keep")]
bm_summary[, rank_high_bm_vw := frank(-avg_bm_vw, ties.method = "average", na.last = "keep")]
bm_summary[, rank_high_bm_median := frank(-avg_bm_median, ties.method = "average", na.last = "keep")]

bm_summary[, rank_high_be_ew := frank(-avg_be_ew, ties.method = "average", na.last = "keep")]
bm_summary[, rank_high_be_vw := frank(-avg_be_vw, ties.method = "average", na.last = "keep")]
bm_summary[, rank_high_be_median := frank(-avg_be_median, ties.method = "average", na.last = "keep")]

cat("\nS14 book-to-market and book-equity rank:\n")
print(
  bm_summary[
    port == NEG_PORT_LABEL,
    .(
      port,
      size,
      op_q,
      inv_q,
      avg_bm_ew,
      avg_bm_vw,
      avg_bm_median,
      rank_high_bm_ew,
      rank_high_bm_vw,
      rank_high_bm_median,
      avg_be_ew,
      avg_be_vw,
      avg_be_median,
      rank_high_be_ew,
      rank_high_be_vw,
      rank_high_be_median,
      avg_n_firms
    )
  ]
)

cat("\nAll 32 portfolios ranked by equal-weighted book-to-market:\n")
print(
  bm_summary[
    order(-avg_bm_ew),
    .(
      port,
      size,
      op_q,
      inv_q,
      avg_bm_ew,
      avg_bm_vw,
      avg_bm_median,
      rank_high_bm_ew,
      avg_n_firms
    )
  ]
)

cat("\nSmall portfolios ranked by equal-weighted book-to-market:\n")
print(
  bm_summary[
    size == "S"
  ][
    order(-avg_bm_ew),
    .(
      port,
      op_q,
      inv_q,
      avg_bm_ew,
      avg_bm_vw,
      avg_bm_median,
      rank_high_bm_ew,
      avg_n_firms
    )
  ]
)

### Industry Composition of S14 ###############################################

s14 <- panel_32[
  port == NEG_PORT_LABEL &
    !is.na(gvkey) &
    !is.na(get(ME_COL)) &
    get(ME_COL) > 0
]

n_months_s14 <- uniqueN(s14$mdate)

### SIC Composition ###########################################################

sic_monthly <- s14[
  !is.na(sic),
  .(
    industry_me = sum(get(ME_COL), na.rm = TRUE),
    n_firms = uniqueN(gvkey),
    n_obs = .N
  ),
  by = .(mdate, sic)
]

sic_monthly[
  ,
  total_me := sum(industry_me, na.rm = TRUE),
  by = mdate
]

sic_monthly[
  ,
  total_firms := sum(n_firms, na.rm = TRUE),
  by = mdate
]

sic_monthly[, me_share := industry_me / total_me]
sic_monthly[, firm_share := n_firms / total_firms]

sic_composition <- sic_monthly[
  ,
  .(
    months_present = uniqueN(mdate),
    avg_me_share = sum(me_share, na.rm = TRUE) / n_months_s14,
    avg_firm_share = sum(firm_share, na.rm = TRUE) / n_months_s14,
    avg_n_firms_when_present = mean_safe(n_firms),
    avg_me_when_present = mean_safe(industry_me)
  ),
  by = sic
][order(-avg_me_share)]

sic_composition[, avg_me_share_pct := 100 * avg_me_share]
sic_composition[, avg_firm_share_pct := 100 * avg_firm_share]
sic_composition[, cumulative_me_share_pct := cumsum(avg_me_share_pct)]
sic_composition[, cumulative_firm_share_pct := cumsum(avg_firm_share_pct)]

cat("\nTop SIC industries in S14 by average monthly market-equity share:\n")
print(sic_composition[1:25])

### SIC2 Composition ##########################################################

sic2_monthly <- s14[
  !is.na(sic2),
  .(
    industry_me = sum(get(ME_COL), na.rm = TRUE),
    n_firms = uniqueN(gvkey),
    n_obs = .N
  ),
  by = .(mdate, sic2)
]

sic2_monthly[
  ,
  total_me := sum(industry_me, na.rm = TRUE),
  by = mdate
]

sic2_monthly[
  ,
  total_firms := sum(n_firms, na.rm = TRUE),
  by = mdate
]

sic2_monthly[, me_share := industry_me / total_me]
sic2_monthly[, firm_share := n_firms / total_firms]

sic2_composition <- sic2_monthly[
  ,
  .(
    months_present = uniqueN(mdate),
    avg_me_share = sum(me_share, na.rm = TRUE) / n_months_s14,
    avg_firm_share = sum(firm_share, na.rm = TRUE) / n_months_s14,
    avg_n_firms_when_present = mean_safe(n_firms),
    avg_me_when_present = mean_safe(industry_me)
  ),
  by = sic2
][order(-avg_me_share)]

sic2_composition[, avg_me_share_pct := 100 * avg_me_share]
sic2_composition[, avg_firm_share_pct := 100 * avg_firm_share]
sic2_composition[, cumulative_me_share_pct := cumsum(avg_me_share_pct)]
sic2_composition[, cumulative_firm_share_pct := cumsum(avg_firm_share_pct)]

cat("\nTop SIC2 industries in S14 by average monthly market-equity share:\n")
print(sic2_composition[1:25])

### Industry Concentration Diagnostics ########################################

sic_hhi_monthly <- sic_monthly[
  ,
  .(
    n_sic = uniqueN(sic),
    hhi_me_sic = sum(me_share^2, na.rm = TRUE),
    hhi_firm_sic = sum(firm_share^2, na.rm = TRUE),
    effective_sic_me = 1 / sum(me_share^2, na.rm = TRUE),
    effective_sic_firm = 1 / sum(firm_share^2, na.rm = TRUE)
  ),
  by = mdate
]

sic2_hhi_monthly <- sic2_monthly[
  ,
  .(
    n_sic2 = uniqueN(sic2),
    hhi_me_sic2 = sum(me_share^2, na.rm = TRUE),
    hhi_firm_sic2 = sum(firm_share^2, na.rm = TRUE),
    effective_sic2_me = 1 / sum(me_share^2, na.rm = TRUE),
    effective_sic2_firm = 1 / sum(firm_share^2, na.rm = TRUE)
  ),
  by = mdate
]

cat("\nS14 SIC concentration summary:\n")
print(
  sic_hhi_monthly[
    ,
    .(
      avg_n_sic = mean_safe(n_sic),
      min_n_sic = min(n_sic, na.rm = TRUE),
      max_n_sic = max(n_sic, na.rm = TRUE),
      avg_hhi_me_sic = mean_safe(hhi_me_sic),
      avg_hhi_firm_sic = mean_safe(hhi_firm_sic),
      avg_effective_sic_me = mean_safe(effective_sic_me),
      avg_effective_sic_firm = mean_safe(effective_sic_firm)
    )
  ]
)

cat("\nS14 SIC2 concentration summary:\n")
print(
  sic2_hhi_monthly[
    ,
    .(
      avg_n_sic2 = mean_safe(n_sic2),
      min_n_sic2 = min(n_sic2, na.rm = TRUE),
      max_n_sic2 = max(n_sic2, na.rm = TRUE),
      avg_hhi_me_sic2 = mean_safe(hhi_me_sic2),
      avg_hhi_firm_sic2 = mean_safe(hhi_firm_sic2),
      avg_effective_sic2_me = mean_safe(effective_sic2_me),
      avg_effective_sic2_firm = mean_safe(effective_sic2_firm)
    )
  ]
)

### Top-Industry Share Summary ################################################

cat("\nShare of S14 accounted for by top SIC industries:\n")
print(
  data.table(
    top_n = c(1, 3, 5, 10, 20),
    me_share_pct = c(
      sic_composition[1, sum(avg_me_share_pct, na.rm = TRUE)],
      sic_composition[1:3, sum(avg_me_share_pct, na.rm = TRUE)],
      sic_composition[1:5, sum(avg_me_share_pct, na.rm = TRUE)],
      sic_composition[1:10, sum(avg_me_share_pct, na.rm = TRUE)],
      sic_composition[1:20, sum(avg_me_share_pct, na.rm = TRUE)]
    ),
    firm_share_pct = c(
      sic_composition[1, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic_composition[1:3, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic_composition[1:5, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic_composition[1:10, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic_composition[1:20, sum(avg_firm_share_pct, na.rm = TRUE)]
    )
  )
)

cat("\nShare of S14 accounted for by top SIC2 industries:\n")
print(
  data.table(
    top_n = c(1, 3, 5, 10, 20),
    me_share_pct = c(
      sic2_composition[1, sum(avg_me_share_pct, na.rm = TRUE)],
      sic2_composition[1:3, sum(avg_me_share_pct, na.rm = TRUE)],
      sic2_composition[1:5, sum(avg_me_share_pct, na.rm = TRUE)],
      sic2_composition[1:10, sum(avg_me_share_pct, na.rm = TRUE)],
      sic2_composition[1:20, sum(avg_me_share_pct, na.rm = TRUE)]
    ),
    firm_share_pct = c(
      sic2_composition[1, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic2_composition[1:3, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic2_composition[1:5, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic2_composition[1:10, sum(avg_firm_share_pct, na.rm = TRUE)],
      sic2_composition[1:20, sum(avg_firm_share_pct, na.rm = TRUE)]
    )
  )
)

### Save Outputs ##############################################################

saveRDS(bm_summary, "S14_OP_INV_2x4x4_bm_rank_summary.rds")
saveRDS(sic_composition, "S14_sic_composition.rds")
saveRDS(sic2_composition, "S14_sic2_composition.rds")
saveRDS(sic_hhi_monthly, "S14_sic_hhi_monthly.rds")
saveRDS(sic2_hhi_monthly, "S14_sic2_hhi_monthly.rds")
