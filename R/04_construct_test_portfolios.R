#### Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
#### Authors: Jacob Korsgaard and Axel Emil Ulvemann
#### Supervisor: Niels Joachim Gormsen
#### Master Thesis 2026

## Timing convention:
## Monthly panel ffyear = t for July t through June t + 1.
## June of calendar year t is used to form portfolios/signals for holding year t.
## Since June itself belongs to ffyear = t - 1 in the monthly panel convention,
## June formation records are mapped using ffyear + 1L when constructing
## holding-period signals and portfolio assignments.

## 4 Test Portfolios ##########################################################
## 4.0 Setup ##################################################################
cat("\014")
rm(list = ls())
graphics.off()

library(data.table)
library(lubridate)
library(zoo)

source(file.path("R", "00_config.R"))

data <- readRDS("crsp_compustat_with_adj_cost.rds")
data[, sic := as.integer(sic)]
if (!"sic2" %in% names(data)) {
  data[, sic2 := fifelse(!is.na(sic), sic %/% 100L, NA_integer_)]
}
setDT(data)

names(data)
uniqueN(data$permno)
uniqueN(data$gvkey)

## 4.1 Formation Pipeline and Timing ##########################################
## Timing convention:
## Monthly panel ffyear = t for July t through June t + 1.
## June of calendar year t is used to form portfolios/signals for holding year t.
## Since June itself belongs to ffyear = t - 1 in the monthly panel convention,
## June formation records are mapped using ffyear + 1L when constructing
## holding-period signals and portfolio assignments.

## 4.1.1 Timing / formation-year conventions ##################################
data[, year := year(mdate)]
data[, month := month(mdate)]
data[, ffyear := fifelse(month >= 7, year, year - 1)]

june_eligible <-
  data[
    month == 6 &
      comp_ok == TRUE &
      has_dec == TRUE &
      has_june == TRUE
  ]

## 4.1.2 Global helper functions ##############################################
ff_20_40_60_80_breaks <- function(x) {
  qs <- as.numeric(quantile(x, probs = c(0.2, 0.4, 0.6, 0.8), type = 1, na.rm = TRUE))
  c(p20 = qs[1], p40 = qs[2], p60 = qs[3], p80 = qs[4])
}

ff_25_50_75_breaks <- function(x) {
  qs <- as.numeric(quantile(x, probs = c(0.25, 0.50, 0.75), type = 1, na.rm = TRUE))
  c(p25 = qs[1], p50 = qs[2], p75 = qs[3])
}

ff_decile_breaks <- function(x) {
  as.numeric(
    quantile(
      x,
      probs = seq(0.1, 0.9, by = 0.1),
      type = 1,
      na.rm = TRUE
    )
  )
}

assign_char4 <- function(x, p25, p50, p75) {
  fcase(
    x <= p25, "1",
    x <= p50, "2",
    x <= p75, "3",
    x > p75, "4",
    default = NA_character_
  )
}

assign_char5 <- function(x, p20, p40, p60, p80, prefix) {
  fcase(
    x <= p20, paste0(prefix, "1"),
    x <= p40, paste0(prefix, "2"),
    x <= p60, paste0(prefix, "3"),
    x <= p80, paste0(prefix, "4"),
    x > p80, paste0(prefix, "5"),
    default = NA_character_
  )
}

wide_to_long_2x4x4 <- function(wide_dt) {
  stopifnot("mdate" %in% names(wide_dt))
  port_cols <- setdiff(names(wide_dt), "mdate")
  long <- melt(
    as.data.table(wide_dt),
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

make_mean_excess_table_2x4x4 <- function(long_dt, rf_dt,
                                         start_date = as.Date("1981-07-01"),
                                         end_date = as.Date("2024-12-01")) {
  x <- long_dt[mdate >= start_date & mdate <= end_date]
  x <- merge(x, rf_dt, by = "mdate", all.x = TRUE)
  x[, ret_excess := ret - RF]

  means <- x[, .(mean_excess = mean(ret_excess, na.rm = TRUE)), by = .(size, c1, c2)]

  grid <- CJ(size = c("S", "B"), c1 = as.character(1:4), c2 = as.character(1:4), unique = TRUE)
  means <- merge(grid, means, by = c("size", "c1", "c2"), all.x = TRUE)

  out <- list()
  for (sz in c("S", "B")) {
    mat <- dcast(means[size == sz], c1 ~ c2, value.var = "mean_excess", drop = FALSE)
    mat_out <- as.matrix(mat[, -1, with = FALSE])
    rownames(mat_out) <- c("Low", "2", "3", "High")
    colnames(mat_out) <- c("Low", "2", "3", "High")
    out[[sz]] <- 100 * mat_out
  }
  out
}

## 4.1.3 FF-style annual signals (OP, INV, Chi) ###############################
## 4.1.3.1 Operating profitability ############################################
op_ff <-
  unique(
    data[, .(gvkey, fyear, revt, cogs, xsga, xint, be_clean)]
  )

op_ff[, op_ff :=
  fifelse(
    !is.na(be_clean) & be_clean > 0 &
      !is.na(revt) & !is.na(cogs),
    (revt
    - cogs
      - fifelse(is.na(xsga), 0, xsga)
      - fifelse(is.na(xint), 0, xint)) / be_clean,
    NA_real_
  )]

op_ff[, fyear := fyear + 1L]

## 4.1.3.2 Investment #########################################################
inv_ff <-
  unique(
    data[, .(gvkey, fyear, at)]
  )

setorder(inv_ff, gvkey, fyear)
inv_ff[, at_lag := shift(at), by = gvkey]
inv_ff[, inv_ff :=
  fifelse(
    !is.na(at_lag) & at_lag > 0,
    (at - at_lag) / at_lag,
    NA_real_
  )]

inv_ff[, fyear := fyear + 1L]

## 4.1.3.3 Adjustment Costs ###################################################
chi_firm_ff <-
  unique(
    data[
      !is.na(gvkey) & !is.na(fyear),
      .(gvkey, fyear, chi_firm, chi_firm_exp)
    ],
    by = c("gvkey", "fyear")
  )
setorder(chi_firm_ff, gvkey, fyear)
chi_firm_ff[, `:=`(
  chi_firm_lag     = shift(chi_firm, 1L),
  chi_firm_exp_lag = shift(chi_firm_exp, 1L)
), by = gvkey]
chi_firm_ff[, fyear := fyear + 1L]

chi_sic_ff <-
  unique(
    data[
      !is.na(sic) & !is.na(fyear),
      .(sic, fyear, chi_sic, chi_sic_exp)
    ],
    by = c("sic", "fyear")
  )
setorder(chi_sic_ff, sic, fyear)
chi_sic_ff[, `:=`(
  chi_sic_lag     = shift(chi_sic, 1L),
  chi_sic_exp_lag = shift(chi_sic_exp, 1L)
), by = sic]
chi_sic_ff[, fyear := fyear + 1L]

chi_sic2_ff <-
  unique(
    data[
      !is.na(sic2) & !is.na(fyear),
      .(sic2, fyear, chi_sic2, chi_sic2_exp)
    ],
    by = c("sic2", "fyear")
  )
setorder(chi_sic2_ff, sic2, fyear)
chi_sic2_ff[, `:=`(
  chi_sic2_lag     = shift(chi_sic2, 1L),
  chi_sic2_exp_lag = shift(chi_sic2_exp, 1L)
), by = sic2]
chi_sic2_ff[, fyear := fyear + 1L]

## 4.1.4 June-t formation table ###############################################
june_form <-
  data[
    month == 6 & has_dec == TRUE & has_june == TRUE,
    .(
      permno,
      gvkey,
      sic,
      sic2,
      ffyear = ffyear + 1L,
      exchcd,
      me_june = me_clean,
      ppe_june = ppegt_lag,
      bm = bm_ff
    )
  ]

june_form <-
  merge(
    june_form,
    op_ff[, .(gvkey, fyear, op_ff)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

june_form <-
  merge(
    june_form,
    inv_ff[, .(gvkey, fyear, inv_ff)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

june_form <-
  merge(
    june_form,
    chi_firm_ff[, .(
      gvkey, fyear,
      chi_firm = chi_firm_lag,
      chi_firm_exp = chi_firm_exp_lag
    )],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE
  )

june_form <-
  merge(
    june_form,
    chi_sic_ff[, .(
      sic, fyear,
      chi_sic = chi_sic_lag,
      chi_sic_exp = chi_sic_exp_lag
    )],
    by.x = c("sic", "ffyear"),
    by.y = c("sic", "fyear"),
    all.x = TRUE
  )

june_form <-
  merge(
    june_form,
    chi_sic2_ff[, .(
      sic2, fyear,
      chi_sic2 = chi_sic2_lag,
      chi_sic2_exp = chi_sic2_exp_lag
    )],
    by.x = c("sic2", "ffyear"),
    by.y = c("sic2", "fyear"),
    all.x = TRUE
  )

setnames(june_form,
  old = c("op_ff", "inv_ff"),
  new = c("op", "inv")
)

setorder(june_form, permno, ffyear)
june_form <- unique(june_form, by = c("permno", "ffyear"))

## 4.1.5 Shared monthly weights ###############################################
setorder(data, permno, mdate)
data[, me_ff_w := shift(me_clean), by = permno]
data[me_ff_w <= 0, me_ff_w := NA_real_]

test_data <-
  data[
    !is.na(retadj) &
      !is.na(me_ff_w) &
      me_ff_w > 0
  ]

out_dir <- INTERIM_DATA_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 4.2 One-dimensional 5 portfolios ##########################################
## Requested families: Size, BM, OP, INV, chi (baseline chi = chi_sic)

## 4.2.1 NYSE breakpoints #####################################################
signals_nyse <- june_form[exchcd == 1]

size1d_bp <-
  signals_nyse[
    !is.na(me_june) & me_june > 0,
    {
      b <- ff_20_40_60_80_breaks(me_june)
      .(size20 = b["p20"], size40 = b["p40"], size60 = b["p60"], size80 = b["p80"])
    },
    by = ffyear
  ]

bm1d_bp <-
  signals_nyse[
    !is.na(bm),
    {
      b <- ff_20_40_60_80_breaks(bm)
      .(bm20 = b["p20"], bm40 = b["p40"], bm60 = b["p60"], bm80 = b["p80"])
    },
    by = ffyear
  ]

op1d_bp <-
  signals_nyse[
    !is.na(op),
    {
      b <- ff_20_40_60_80_breaks(op)
      .(op20 = b["p20"], op40 = b["p40"], op60 = b["p60"], op80 = b["p80"])
    },
    by = ffyear
  ]

inv1d_bp <-
  signals_nyse[
    !is.na(inv),
    {
      b <- ff_20_40_60_80_breaks(inv)
      .(inv20 = b["p20"], inv40 = b["p40"], inv60 = b["p60"], inv80 = b["p80"])
    },
    by = ffyear
  ]

chi1d_bp <-
  signals_nyse[
    !is.na(chi_sic),
    {
      b <- ff_20_40_60_80_breaks(chi_sic)
      .(chi20 = b["p20"], chi40 = b["p40"], chi60 = b["p60"], chi80 = b["p80"])
    },
    by = ffyear
  ]

june_form <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
    list(june_form, size1d_bp, bm1d_bp, op1d_bp, inv1d_bp, chi1d_bp)
  )

## 4.2.2 Assign labels ########################################################
june_form[, size1d := assign_char5(me_june, size20, size40, size60, size80, "S")]
june_form[, bm1d := assign_char5(bm, bm20, bm40, bm60, bm80, "B")]
june_form[, op1d := assign_char5(op, op20, op40, op60, op80, "O")]
june_form[, inv1d := assign_char5(inv, inv20, inv40, inv60, inv80, "I")]
june_form[, chi1d := assign_char5(chi_sic, chi20, chi40, chi60, chi80, "C")]

## 4.2.3 Merge labels into monthly panel ######################################
data <-
  merge(
    data,
    june_form[, .(permno, ffyear, size1d, bm1d, op1d, inv1d, chi1d)],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

test_data <-
  data[
    !is.na(retadj) &
      !is.na(me_ff_w) &
      me_ff_w > 0
  ]

## 4.2.4 Builder ##############################################################
build_1d_5port <- function(dt, label_var, file_stub) {
  port_dt <-
    dt[
      !is.na(get(label_var)),
      .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
      by = .(mdate, bucket = get(label_var))
    ]

  port_dt[, port := bucket]

  wide <-
    dcast(
      port_dt[, .(mdate, port, vwret)],
      mdate ~ port,
      value.var = "vwret"
    )

  setcolorder(wide, c("mdate", sort(setdiff(names(wide), "mdate"))))

  saveRDS(
    wide,
    file = file.path(out_dir, paste0("my_", file_stub, "_1D_wide.rds"))
  )

  list(long = port_dt, wide = wide)
}

## 4.2.5 Run one-dimensional 5 sorts ##########################################
res_SIZE_1D <- build_1d_5port(test_data, "size1d", "SIZE")
res_BM1D_1D <- build_1d_5port(test_data, "bm1d", "BM1D")
res_OP1D_1D <- build_1d_5port(test_data, "op1d", "OP1D")
res_INV1D_1D <- build_1d_5port(test_data, "inv1d", "INV1D")
res_CHI1D_1D <- build_1d_5port(test_data, "chi1d", "CHI1D")

cat("\nSaved 1D 5-portfolio panels:\n")
cat("my_SIZE_1D_wide.rds\n")
cat("my_BM1D_1D_wide.rds\n")
cat("my_OP1D_1D_wide.rds\n")
cat("my_INV1D_1D_wide.rds\n")
cat("my_CHI1D_1D_wide.rds\n")

## 4.3 5x5 Portfolios #########################################################
## Existing families: Size x BM, Size x OP, Size x INV, Size x chi variants
## New requested families: BM x chi, OP x chi, INV x chi (baseline chi = chi_sic)

## 4.3.1 NYSE universe and 5-bin breakpoints ##################################
signals_nyse <- june_form[exchcd == 1]

size_bp <-
  signals_nyse[
    !is.na(me_june) & me_june > 0,
    {
      b <- ff_20_40_60_80_breaks(me_june)
      .(me20 = b["p20"], me40 = b["p40"], me60 = b["p60"], me80 = b["p80"])
    },
    by = ffyear
  ]


chi_vars <- c(
  "chi_firm", "chi_firm_exp",
  "chi_sic",  "chi_sic_exp",
  "chi_sic2", "chi_sic2_exp"
)

chi_bp_list <-
  lapply(chi_vars, function(v) {
    signals_nyse[
      !is.na(get(v)),
      {
        b <- ff_20_40_60_80_breaks(get(v))
        setNames(
          as.list(b),
          paste0(v, c("20", "40", "60", "80"))
        )
      },
      by = ffyear
    ]
  })

chi_bp <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all = TRUE),
    chi_bp_list
  )

june_form <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
    list(june_form, size_bp, chi_bp)
  )

## 4.3.2 Existing 5x5 labels ##################################################
june_form[, size5 :=
  fcase(
    me_june <= me20, "S1",
    me_june <= me40, "S2",
    me_june <= me60, "S3",
    me_june <= me80, "S4",
    me_june > me80, "S5",
    default = NA_character_
  )]

june_form[, value5 :=
  fcase(
    bm <= bm20, "V1",
    bm <= bm40, "V2",
    bm <= bm60, "V3",
    bm <= bm80, "V4",
    bm > bm80, "V5",
    default = NA_character_
  )]

june_form[, prof5 :=
  fcase(
    op <= op20, "P1",
    op <= op40, "P2",
    op <= op60, "P3",
    op <= op80, "P4",
    op > op80, "P5",
    default = NA_character_
  )]

june_form[, invest5 :=
  fcase(
    inv <= inv20, "I1",
    inv <= inv40, "I2",
    inv <= inv60, "I3",
    inv <= inv80, "I4",
    inv > inv80, "I5",
    default = NA_character_
  )]

for (v in chi_vars) {
  q20 <- paste0(v, "20")
  q40 <- paste0(v, "40")
  q60 <- paste0(v, "60")
  q80 <- paste0(v, "80")
  lab <- paste0(v, "5")

  june_form[, (lab) :=
    fcase(
      get(v) <= get(q20), "C1",
      get(v) <= get(q40), "C2",
      get(v) <= get(q60), "C3",
      get(v) <= get(q80), "C4",
      get(v) > get(q80), "C5",
      default = NA_character_
    )]
}

## 4.3.3 New requested 5x5 labels #############################################
## Baseline chi interaction uses chi_sic only
june_form[, bm5_chi := assign_char5(bm, bm20, bm40, bm60, bm80, "B")]
june_form[, op5_chi := assign_char5(op, op20, op40, op60, op80, "O")]
june_form[, inv5_chi := assign_char5(inv, inv20, inv40, inv60, inv80, "I")]
june_form[, chi5 := assign_char5(chi_sic, chi_sic20, chi_sic40, chi_sic60, chi_sic80, "C")]

label_vars <- c(
  "size5", "value5", "prof5", "invest5",
  paste0(chi_vars, "5"),
  "bm5_chi", "op5_chi", "inv5_chi", "chi5",
  "me_june", "ppe_june"
)

data <-
  merge(
    data,
    june_form[, c("permno", "ffyear", label_vars), with = FALSE],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

setorder(data, permno, mdate)
data[, me_ff_w := shift(me_clean), by = permno]
data[me_ff_w <= 0, me_ff_w := NA_real_]

test_data <-
  data[
    !is.na(retadj) &
      !is.na(me_ff_w) &
      me_ff_w > 0
  ]

## 4.3.4 Existing 5x5 portfolio construction ##################################
port_5x5 <-
  test_data[
    !is.na(size5) & !is.na(value5),
    .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
    by = .(mdate, size5, value5)
  ]
port_5x5[, port := paste0(size5, value5)]

port_5x5_wide <-
  dcast(
    port_5x5,
    mdate ~ port,
    value.var = "vwret"
  )

port_5x5_op <-
  test_data[
    !is.na(size5) & !is.na(prof5),
    .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
    by = .(mdate, size5, prof5)
  ]
port_5x5_op[, port := paste0(size5, prof5)]

port_5x5_op_wide <-
  dcast(
    port_5x5_op,
    mdate ~ port,
    value.var = "vwret"
  )

port_5x5_inv <-
  test_data[
    !is.na(size5) & !is.na(invest5),
    .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
    by = .(mdate, size5, invest5)
  ]
port_5x5_inv[, port := paste0(size5, invest5)]

port_5x5_inv_wide <-
  dcast(
    port_5x5_inv,
    mdate ~ port,
    value.var = "vwret"
  )

size_chi_ports <- list()
size_chi_wide <- list()

for (v in chi_vars) {
  chi_lab <- paste0(v, "5")

  port_dt <-
    test_data[
      !is.na(size5) & !is.na(get(chi_lab)),
      .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
      by = .(mdate, size5, chi = get(chi_lab))
    ]

  port_dt[, port := paste0(size5, chi)]
  size_chi_ports[[v]] <- port_dt

  size_chi_wide[[v]] <-
    dcast(
      port_dt,
      mdate ~ port,
      value.var = "vwret"
    )
}

## 4.3.5 New requested 5x5 families: BM x chi, OP x chi, INV x chi ###########
build_5x5_family <- function(dt, row_var, col_var, file_stub) {
  port_dt <-
    dt[
      !is.na(get(row_var)) & !is.na(get(col_var)),
      .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
      by = .(mdate, row = get(row_var), col = get(col_var))
    ]

  port_dt[, port := paste0(row, col)]

  wide <-
    dcast(
      port_dt[, .(mdate, port, vwret)],
      mdate ~ port,
      value.var = "vwret"
    )

  setcolorder(wide, c("mdate", sort(setdiff(names(wide), "mdate"))))

  saveRDS(
    wide,
    file = file.path(out_dir, paste0("my_", file_stub, "_5x5_wide.rds"))
  )

  list(long = port_dt, wide = wide)
}

res_BM_CHI_SIC_5x5 <- build_5x5_family(test_data, "bm5_chi", "chi5", "BM_CHI_SIC")
res_OP_CHI_SIC_5x5 <- build_5x5_family(test_data, "op5_chi", "chi5", "OP_CHI_SIC")
res_INV_CHI_SIC_5x5 <- build_5x5_family(test_data, "inv5_chi", "chi5", "INV_CHI_SIC")

## 4.3.6 Mean/SD/SR for existing 5x5 ##########################################
START_DATE <- as.Date("1981-07-01")
END_DATE <- as.Date("2024-12-01")

filter_window <- function(dt) dt[mdate >= START_DATE & mdate <= END_DATE]

make_5x5_mean_matrix <- function(port_dt,
                                 size_var = "size5",
                                 char_var,
                                 ret_var = "vwret") {
  means <-
    port_dt[
      ,
      .(mean_ret = mean(get(ret_var), na.rm = TRUE)),
      by = c(size_var, char_var)
    ]

  means[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  means[, (char_var) := factor(get(char_var),
    levels = sort(unique(get(char_var)))
  )]

  mat <-
    dcast(means,
      as.formula(paste(size_var, "~", char_var)),
      value.var = "mean_ret"
    )

  mat_out <- as.matrix(mat[, -1])
  rownames(mat_out) <- mat[[size_var]]

  100 * mat_out
}

make_5x5_sd_matrix <- function(port_dt,
                               size_var = "size5",
                               char_var,
                               ret_var = "vwret_excess") {
  sds <-
    port_dt[
      ,
      .(sd_ret = sd(get(ret_var), na.rm = TRUE)),
      by = c(size_var, char_var)
    ]

  sds[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  sds[, (char_var) := factor(get(char_var),
    levels = sort(unique(get(char_var)))
  )]

  mat <-
    dcast(sds,
      as.formula(paste(size_var, "~", char_var)),
      value.var = "sd_ret"
    )

  mat_out <- as.matrix(mat[, -1])
  rownames(mat_out) <- mat[[size_var]]

  100 * mat_out
}

make_5x5_sr_matrix <- function(port_dt,
                               size_var = "size5",
                               char_var,
                               ret_var = "vwret_excess") {
  stats <-
    port_dt[
      ,
      .(
        mean_ret = mean(get(ret_var), na.rm = TRUE),
        sd_ret   = sd(get(ret_var), na.rm = TRUE)
      ),
      by = c(size_var, char_var)
    ]

  stats[, sr := mean_ret / sd_ret]

  stats[, (size_var) := factor(get(size_var), levels = paste0("S", 1:5))]
  stats[, (char_var) := factor(get(char_var),
    levels = sort(unique(get(char_var)))
  )]

  mat <-
    dcast(stats,
      as.formula(paste(size_var, "~", char_var)),
      value.var = "sr"
    )

  mat_out <- as.matrix(mat[, -1])
  rownames(mat_out) <- mat[[size_var]]

  mat_out
}

mean_5x5 <- list()
mean_5x5[["bm"]] <- make_5x5_mean_matrix(port_5x5, "size5", "value5")
mean_5x5[["op"]] <- make_5x5_mean_matrix(port_5x5_op, "size5", "prof5")
mean_5x5[["inv"]] <- make_5x5_mean_matrix(port_5x5_inv, "size5", "invest5")

for (v in names(size_chi_ports)) {
  mean_5x5[[v]] <- make_5x5_mean_matrix(size_chi_ports[[v]], "size5", "chi")
}

rf_monthly <- unique(data[, .(mdate, RF)])

make_excess <- function(port_dt) {
  tmp <- merge(port_dt, rf_monthly, by = "mdate", all.x = TRUE)
  tmp[, vwret_excess := vwret - RF]
  tmp
}

mean_excess_5x5 <- list()
sd_5x5 <- list()
sr_5x5 <- list()

build_excess_stats <- function(port_dt, char_var) {
  ex <- make_excess(filter_window(port_dt))
  list(
    mean = make_5x5_mean_matrix(ex, "size5", char_var, "vwret_excess"),
    sd   = make_5x5_sd_matrix(ex, "size5", char_var, "vwret_excess"),
    sr   = make_5x5_sr_matrix(ex, "size5", char_var, "vwret_excess")
  )
}

tmp <- build_excess_stats(port_5x5, "value5")
mean_excess_5x5[["bm"]] <- tmp$mean
sd_5x5[["bm"]] <- tmp$sd
sr_5x5[["bm"]] <- tmp$sr

tmp <- build_excess_stats(port_5x5_op, "prof5")
mean_excess_5x5[["op"]] <- tmp$mean
sd_5x5[["op"]] <- tmp$sd
sr_5x5[["op"]] <- tmp$sr

tmp <- build_excess_stats(port_5x5_inv, "invest5")
mean_excess_5x5[["inv"]] <- tmp$mean
sd_5x5[["inv"]] <- tmp$sd
sr_5x5[["inv"]] <- tmp$sr

for (v in names(size_chi_ports)) {
  tmp <- build_excess_stats(size_chi_ports[[v]], "chi")
  mean_excess_5x5[[v]] <- tmp$mean
  sd_5x5[[v]] <- tmp$sd
  sr_5x5[[v]] <- tmp$sr
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

display_all_5x5(mean_5x5, "MONTHLY MEAN 5x5 RETURNS (%)")
display_all_5x5(mean_excess_5x5, "MONTHLY MEAN 5x5 EXCESS RETURNS")
display_all_5x5(sd_5x5, "MONTHLY STD DEV 5x5 EXCESS RETURNS (%)")
display_all_5x5(sr_5x5, "SHARPE RATIOS 5x5 (monthly)")

## 4.3.7 Save existing and new 5x5 panels #####################################
saveRDS(port_5x5_wide, file = file.path(out_dir, "my_Size_BM_5x5_wide.rds"))
saveRDS(port_5x5_op_wide, file = file.path(out_dir, "my_Size_OP_5x5_wide.rds"))
saveRDS(port_5x5_inv_wide, file = file.path(out_dir, "my_Size_INV_5x5_wide.rds"))

if (exists("size_chi_wide") && is.list(size_chi_wide) && length(size_chi_wide) > 0) {
  for (v in names(size_chi_wide)) {
    saveRDS(
      size_chi_wide[[v]],
      file = file.path(out_dir, paste0("my_Size_", v, "_5x5_wide.rds"))
    )
  }
  cat("Also saved chi 5x5 wide panels as RDS (my_Size_<chi>_5x5_wide.rds)\n\n")
}

cat("Saved new requested 5x5 chi interaction panels:\n")
cat("my_BM_CHI_SIC_5x5_wide.rds\n")
cat("my_OP_CHI_SIC_5x5_wide.rds\n")
cat("my_INV_CHI_SIC_5x5_wide.rds\n\n")

## 4.4 2x4x4 Portfolios #######################################################
## 4.4.0 Output directory #####################################################
out_dir <- INTERIM_DATA_DIR
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 4.4.1 Helper for 2x4x4 #####################################################
wide_to_long_2x4x4 <- wide_to_long_2x4x4

## 4.4.2 Core builder #########################################################
build_2x4x4_family <- function(june_form, data, var1, var2, name_tag,
                               start_date = as.Date("1981-07-01"),
                               end_date = as.Date("2024-12-01"),
                               save_rds = TRUE) {
  stopifnot(june_form[, uniqueN(paste(permno, ffyear))] == nrow(june_form))

  rf_monthly <- unique(data[!is.na(RF), .(mdate = as.Date(mdate), RF)])

  nyse <- june_form[exchcd == 1 & !is.na(me_june) & me_june > 0]

  size_bp <- nyse[, .(me50 = median(me_june, na.rm = TRUE)), by = ffyear]

  nyse <- merge(nyse, size_bp, by = "ffyear", all.x = TRUE)
  nyse[, size2 := fifelse(me_june <= me50, "S", "B")]

  bp_v1 <-
    nyse[
      !is.na(get(var1)),
      {
        b <- ff_25_50_75_breaks(get(var1))
        .(v1_25 = b["p25"], v1_50 = b["p50"], v1_75 = b["p75"])
      },
      by = .(ffyear, size2)
    ]

  bp_v2 <-
    nyse[
      !is.na(get(var2)),
      {
        b <- ff_25_50_75_breaks(get(var2))
        .(v2_25 = b["p25"], v2_50 = b["p50"], v2_75 = b["p75"])
      },
      by = .(ffyear, size2)
    ]

  jf <- merge(june_form, size_bp, by = "ffyear", all.x = TRUE)
  jf[, size2 := fifelse(me_june <= me50, "S", "B")]

  jf <- merge(jf, bp_v1, by = c("ffyear", "size2"), all.x = TRUE)
  jf <- merge(jf, bp_v2, by = c("ffyear", "size2"), all.x = TRUE)

  jf[, c1 := assign_char4(get(var1), v1_25, v1_50, v1_75)]
  jf[, c2 := assign_char4(get(var2), v2_25, v2_50, v2_75)]

  port_col <- paste0("port_", name_tag, "_2x4x4")
  jf[
    !is.na(size2) & !is.na(c1) & !is.na(c2),
    (port_col) := paste0(size2, c1, c2)
  ]

  stopifnot(jf[, uniqueN(paste(permno, ffyear))] == nrow(jf))

  labels <-
    jf[
      !is.na(get(port_col)),
      .(permno, ffyear, port = get(port_col))
    ]

  stopifnot(labels[, uniqueN(paste(permno, ffyear))] == nrow(labels))

  labels_named <-
    jf[
      ,
      .(permno, ffyear, label = get(port_col))
    ]
  setnames(labels_named, "label", port_col)

  panel <- merge(data, labels, by = c("permno", "ffyear"), all.x = TRUE)

  ports <-
    panel[
      !is.na(port) &
        !is.na(me_ff_w) & me_ff_w > 0 &
        !is.na(retadj),
      .(vwret = sum(me_ff_w * retadj) / sum(me_ff_w)),
      by = .(mdate, port)
    ]

  ports_wide <- dcast(ports, mdate ~ port, value.var = "vwret")

  if (save_rds) {
    rds_name <- paste0("my_Size_", toupper(name_tag), "_2x4x4_wide.rds")
    saveRDS(ports_wide, file = file.path(out_dir, rds_name))

    label_name <- paste0("labels_", toupper(name_tag), "_2x4x4.rds")
    saveRDS(labels_named, file = file.path(out_dir, label_name))
  }

  long <- wide_to_long_2x4x4(ports_wide)
  mean_excess <- make_mean_excess_table_2x4x4(long, rf_monthly, start_date, end_date)

  list(
    wide = ports_wide,
    mean_excess = mean_excess,
    port_col = port_col,
    labels_named = labels_named
  )
}

## 4.4.3 Run 2x4x4 families ###################################################
portfolio_families <- list(
  list(var1 = "bm", var2 = "op", tag = "BM_OP"),
  list(var1 = "bm", var2 = "inv", tag = "BM_INV"),
  list(var1 = "op", var2 = "inv", tag = "OP_INV"),
  list(var1 = "op", var2 = "chi_firm", tag = "OP_CHI_FIRM"),
  list(var1 = "op", var2 = "chi_sic", tag = "OP_CHI_SIC"),
  list(var1 = "op", var2 = "chi_sic2", tag = "OP_CHI_SIC2"),
  list(var1 = "bm", var2 = "chi_firm", tag = "BM_CHI_FIRM"),
  list(var1 = "bm", var2 = "chi_sic", tag = "BM_CHI_SIC"),
  list(var1 = "bm", var2 = "chi_sic2", tag = "BM_CHI_SIC2"),
  list(var1 = "inv", var2 = "chi_firm", tag = "INV_CHI_FIRM"),
  list(var1 = "inv", var2 = "chi_sic", tag = "INV_CHI_SIC"),
  list(var1 = "inv", var2 = "chi_sic2", tag = "INV_CHI_SIC2")
)

results_2x4x4 <- list()

for (pf in portfolio_families) {
  cat("\n\n=============================================================\n")
  cat("AVERAGE MONTHLY PERCENT EXCESS RETURNS\n")
  cat("2x4x4:", pf$tag, "\n")
  cat("Sample: 1981-07-01 to 2024-12-01\n")
  cat("=============================================================\n\n")

  res <-
    build_2x4x4_family(
      june_form  = june_form,
      data       = data,
      var1       = pf$var1,
      var2       = pf$var2,
      name_tag   = pf$tag,
      start_date = as.Date("1981-07-01"),
      end_date   = as.Date("2024-12-01"),
      save_rds   = TRUE
    )

  results_2x4x4[[pf$tag]] <- res

  cat("Small\n")
  print(round(res$mean_excess$S, 2))

  cat("\nBig\n")
  print(round(res$mean_excess$B, 2))
}

## Merge all 2x4x4 labels into main data ######################################
labels_2x4x4_list <-
  lapply(results_2x4x4, function(x) x$labels_named)

labels_2x4x4_merged <-
  Reduce(
    function(x, y) merge(x, y, by = c("permno", "ffyear"), all = TRUE),
    labels_2x4x4_list
  )

data <-
  merge(
    data,
    labels_2x4x4_merged,
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

cat("\nMerged 2x4x4 family labels into main data.\n")
print(grep("^port_.*_2x4x4$", names(data), value = TRUE))

## 4.5 Decile Portfolios for Chi ##############################################
chi_vars <- c(
  "chi_firm", "chi_firm_exp",
  "chi_sic",  "chi_sic_exp",
  "chi_sic2", "chi_sic2_exp"
)

signals_nyse <- june_form[exchcd == 1]

chi10_bp_list <-
  lapply(chi_vars, function(v) {
    signals_nyse[
      !is.na(get(v)),
      {
        b <- ff_decile_breaks(get(v))
        as.list(setNames(b, paste0(v, "_p", seq(10, 90, by = 10))))
      },
      by = ffyear
    ]
  })

chi10_bp <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all = TRUE),
    chi10_bp_list
  )

june_form <-
  merge(
    june_form,
    chi10_bp,
    by = "ffyear",
    all.x = TRUE
  )

for (v in chi_vars) {
  june_form[, paste0(v, "_dec10") :=
    fcase(
      get(v) <= get(paste0(v, "_p10")), "D1",
      get(v) <= get(paste0(v, "_p20")), "D2",
      get(v) <= get(paste0(v, "_p30")), "D3",
      get(v) <= get(paste0(v, "_p40")), "D4",
      get(v) <= get(paste0(v, "_p50")), "D5",
      get(v) <= get(paste0(v, "_p60")), "D6",
      get(v) <= get(paste0(v, "_p70")), "D7",
      get(v) <= get(paste0(v, "_p80")), "D8",
      get(v) <= get(paste0(v, "_p90")), "D9",
      get(v) > get(paste0(v, "_p90")), "D10",
      default = NA_character_
    )]
}

dec10_vars <- paste0(chi_vars, "_dec10")

data <-
  merge(
    data,
    june_form[, c("permno", "ffyear", dec10_vars), with = FALSE],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

cat("\nChecking final labelled variables before save:\n")
print(grep("^port_.*_2x4x4$|_dec10$|5$|1d$", names(data), value = TRUE, ignore.case = TRUE))

saveRDS(data, "crsp_compustat_with_portfolio_labels.rds")
cat("\nSaved: crsp_compustat_with_portfolio_labels.rds\n")

chi10_results <- list()

for (v in chi_vars) {
  dec_var <- paste0(v, "_dec10")

  port_dt <-
    data[
      !is.na(get(dec_var)) &
        !is.na(me_ff_w) & me_ff_w > 0 &
        !is.na(retadj),
      .(
        vwret   = sum(me_ff_w * retadj) / sum(me_ff_w),
        n_firms = uniqueN(permno)
      ),
      by = .(mdate, dec10 = get(dec_var))
    ]

  setorder(port_dt, mdate, dec10)

  port_wide <-
    dcast(
      port_dt,
      mdate ~ dec10,
      value.var = "vwret"
    )

  saveRDS(
    port_wide,
    file = file.path(out_dir, paste0("my_", v, "_DEC10_wide.rds"))
  )

  chi10_results[[v]] <- list(
    long = port_dt,
    wide = port_wide
  )
}

for (v in names(chi10_results)) {
  cat("\n============================================================\n")
  cat("CHI DECILES (1×10, JUNE-FORMED)\n")
  cat("Sorting variable:", v, "\n")
  cat("============================================================\n\n")

  dt <- chi10_results[[v]]$long

  means <-
    dt[
      ,
      .(mean_ret = mean(vwret, na.rm = TRUE)),
      by = dec10
    ][order(dec10)]

  means[, mean_ret_pct := 100 * mean_ret]

  print(means[, .(dec10, mean_ret_pct)])
}


june_form[, .N, by = .(permno, ffyear)][N > 1]


data[, .N, by = .(permno, mdate)][N > 1]
names(june_form)[grepl("^(bm|op|inv|me|chi)", names(june_form))]
list.files(out_dir, pattern = "my_.*(1D|5x5|DEC10).*\\.rds$")

file.exists(file.path(out_dir, "my_Size_OP_5x5_wide.rds"))
exists("port_5x5_op_wide")
dim(port_5x5_op_wide)
names(port_5x5_op_wide)


## 4.6 Long-Short Chi Portfolios ##############################################

rf_monthly <- unique(data[!is.na(RF), .(mdate = as.Date(mdate), RF)])

## 4.6.1 Size x Chi: 5 long-short series (one per size quintile) ##############
## Long C5 (high chi = slow), Short C1 (low chi = fast)

size_chi_ls <- list()

for (sz in paste0("S", 1:5)) {
  dt <- size_chi_ports[["chi_sic"]][size5 == sz]

  long_leg <- dt[chi == "C5", .(mdate, ret_long = vwret)]
  short_leg <- dt[chi == "C1", .(mdate, ret_short = vwret)]

  ls <- merge(long_leg, short_leg, by = "mdate", all = TRUE)
  ls[, ls_ret := ret_long - ret_short]
  ls[, size := sz]

  size_chi_ls[[sz]] <- ls
}

size_chi_ls_dt <- rbindlist(size_chi_ls)

size_chi_ls_wide <- dcast(
  size_chi_ls_dt[, .(mdate, size, ls_ret)],
  mdate ~ size,
  value.var = "ls_ret"
)
setnames(size_chi_ls_wide,
  old = paste0("S", 1:5),
  new = paste0("LS_chi_", paste0("S", 1:5))
)

saveRDS(size_chi_ls_wide,
  file = file.path(out_dir, "my_Size_CHI_LS_wide.rds")
)

## 4.6.2 2x4x4 long-short series: 8 series each for OP, BM, INV ##############
## For each (size2, c1) cell: long C4 (high chi), short C1 (low chi)

build_2x4x4_ls <- function(ports_wide, name_tag) {
  ## Convert to long with size, c1, c2 columns
  long <- wide_to_long_2x4x4(ports_wide)

  ## c2 is the chi dimension; long C4, short C1
  long_leg <- long[c2 == "4", .(mdate, size, c1, ret_long = ret)]
  short_leg <- long[c2 == "1", .(mdate, size, c1, ret_short = ret)]

  ls <- merge(long_leg, short_leg, by = c("mdate", "size", "c1"), all = TRUE)
  ls[, ls_ret := ret_long - ret_short]
  ls[, port_id := paste0("LS_", name_tag, "_", size, c1)]

  wide <- dcast(
    ls[, .(mdate, port_id, ls_ret)],
    mdate ~ port_id,
    value.var = "ls_ret"
  )

  saveRDS(wide,
    file = file.path(out_dir, paste0("my_", name_tag, "_LS_wide.rds"))
  )

  list(long = ls, wide = wide)
}

ls_OP_CHI_SIC <- build_2x4x4_ls(results_2x4x4[["OP_CHI_SIC"]]$wide, "OP_CHI_SIC")
ls_BM_CHI_SIC <- build_2x4x4_ls(results_2x4x4[["BM_CHI_SIC"]]$wide, "BM_CHI_SIC")
ls_INV_CHI_SIC <- build_2x4x4_ls(results_2x4x4[["INV_CHI_SIC"]]$wide, "INV_CHI_SIC")

## 4.6.3 Summary: mean excess returns for all long-short series ################

summarise_ls <- function(ls_dt, tag,
                         start_date = as.Date("1981-07-01"),
                         end_date = as.Date("2024-12-01")) {
  dt <- ls_dt[mdate >= start_date & mdate <= end_date]
  dt <- merge(dt, rf_monthly, by = "mdate", all.x = TRUE)
  dt[, ls_excess := ls_ret - RF]

  stats <- dt[,
    .(
      mean_excess = round(100 * mean(ls_excess, na.rm = TRUE), 3),
      sd = round(100 * sd(ls_excess, na.rm = TRUE), 3),
      sr = round(mean(ls_excess, na.rm = TRUE) /
        sd(ls_excess, na.rm = TRUE), 3)
    ),
    by = port_id
  ]

  cat("\n=============================================\n")
  cat("LONG-SHORT SUMMARY:", tag, "\n")
  cat("Sample: 1981-07-01 to 2024-12-01\n")
  cat("=============================================\n")
  print(stats[order(port_id)])
}

## Size x Chi summary
size_chi_ls_summary <- size_chi_ls_dt[mdate >= as.Date("1981-07-01") &
  mdate <= as.Date("2024-12-01")]
size_chi_ls_summary <- merge(size_chi_ls_summary, rf_monthly, by = "mdate", all.x = TRUE)
size_chi_ls_summary[, ls_excess := ls_ret - RF]

cat("\n=============================================\n")
cat("LONG-SHORT SUMMARY: Size x Chi\n")
cat("Sample: 1981-07-01 to 2024-12-01\n")
cat("=============================================\n")
print(
  size_chi_ls_summary[,
    .(
      mean_excess = round(100 * mean(ls_excess, na.rm = TRUE), 3),
      sd = round(100 * sd(ls_excess, na.rm = TRUE), 3),
      sr = round(mean(ls_excess, na.rm = TRUE) /
        sd(ls_excess, na.rm = TRUE), 3)
    ),
    by = size
  ][order(size)]
)

summarise_ls(ls_OP_CHI_SIC$long, "Size x OP x Chi")
summarise_ls(ls_BM_CHI_SIC$long, "Size x BM x Chi")
summarise_ls(ls_INV_CHI_SIC$long, "Size x INV x Chi")

## 4.7 Additional Long-Short Portfolios #######################################

## Conventional directions:
## Size: S1 - S5 (small minus big)
## BM:   High - Low (V5 - V1)
## OP:   High - Low (P5 - P1)
## INV:  Low - High (I1 - I5)
## Chi:  High - Low (C5 - C1) -- already done in 4.6

## 4.7.1 Size long-short within each Chi quintile (from 5x5 Size x Chi) #######

size_ls_within_chi <- list()

for (chi_q in paste0("C", 1:5)) {
  dt <- size_chi_ports[["chi_sic"]][chi == chi_q]

  long_leg <- dt[size5 == "S1", .(mdate, ret_long = vwret)]
  short_leg <- dt[size5 == "S5", .(mdate, ret_short = vwret)]

  ls <- merge(long_leg, short_leg, by = "mdate", all = TRUE)
  ls[, ls_ret := ret_long - ret_short]
  ls[, port_id := paste0("LS_Size_", chi_q)]

  size_ls_within_chi[[chi_q]] <- ls
}

size_ls_within_chi_dt <- rbindlist(size_ls_within_chi)

size_ls_within_chi_wide <- dcast(
  size_ls_within_chi_dt[, .(mdate, port_id, ls_ret)],
  mdate ~ port_id,
  value.var = "ls_ret"
)

saveRDS(size_ls_within_chi_wide,
  file = file.path(out_dir, "my_Size_LS_within_CHI_wide.rds")
)

## Summary
size_ls_within_chi_summary <-
  size_ls_within_chi_dt[mdate >= as.Date("1981-07-01") &
    mdate <= as.Date("2024-12-01")]
size_ls_within_chi_summary <-
  merge(size_ls_within_chi_summary, rf_monthly, by = "mdate", all.x = TRUE)
size_ls_within_chi_summary[, ls_excess := ls_ret - RF]

cat("\n=============================================\n")
cat("LONG-SHORT SUMMARY: Size (S1-S5) within Chi quintile\n")
cat("Sample: 1981-07-01 to 2024-12-01\n")
cat("=============================================\n")
print(
  size_ls_within_chi_summary[,
    .(
      mean_excess = round(100 * mean(ls_excess, na.rm = TRUE), 3),
      sd = round(100 * sd(ls_excess, na.rm = TRUE), 3),
      sr = round(mean(ls_excess, na.rm = TRUE) /
        sd(ls_excess, na.rm = TRUE), 3)
    ),
    by = port_id
  ][order(port_id)]
)

## 4.7.2 Characteristic long-short within each Chi quartile (2x4x4) ###########
## For each (size2, c2) cell: long-short on c1 using conventional direction
## BM:  c1 == "4" minus c1 == "1"
## OP:  c1 == "4" minus c1 == "1"
## INV: c1 == "1" minus c1 == "4"

build_2x4x4_char_ls <- function(ports_wide, name_tag, long_c1, short_c1) {
  long_dt <- wide_to_long_2x4x4(ports_wide)

  long_leg <- long_dt[c1 == long_c1, .(mdate, size, c2, ret_long = ret)]
  short_leg <- long_dt[c1 == short_c1, .(mdate, size, c2, ret_short = ret)]

  ls <- merge(long_leg, short_leg, by = c("mdate", "size", "c2"), all = TRUE)
  ls[, ls_ret := ret_long - ret_short]
  ls[, port_id := paste0("LS_", name_tag, "_", size, "_Chi", c2)]

  wide <- dcast(
    ls[, .(mdate, port_id, ls_ret)],
    mdate ~ port_id,
    value.var = "ls_ret"
  )

  saveRDS(wide,
    file = file.path(out_dir, paste0("my_", name_tag, "_char_LS_wide.rds"))
  )

  list(long = ls, wide = wide)
}

ls_char_OP_CHI_SIC <-
  build_2x4x4_char_ls(
    results_2x4x4[["OP_CHI_SIC"]]$wide,
    name_tag  = "OP_CHI_SIC",
    long_c1   = "4",
    short_c1  = "1"
  )

ls_char_BM_CHI_SIC <-
  build_2x4x4_char_ls(
    results_2x4x4[["BM_CHI_SIC"]]$wide,
    name_tag  = "BM_CHI_SIC",
    long_c1   = "4",
    short_c1  = "1"
  )

ls_char_INV_CHI_SIC <-
  build_2x4x4_char_ls(
    results_2x4x4[["INV_CHI_SIC"]]$wide,
    name_tag  = "INV_CHI_SIC",
    long_c1   = "1", ## INV: low minus high
    short_c1  = "4"
  )

## Summaries
summarise_ls(ls_char_OP_CHI_SIC$long, "OP (H-L) within Chi quartile")
summarise_ls(ls_char_BM_CHI_SIC$long, "BM (H-L) within Chi quartile")
summarise_ls(ls_char_INV_CHI_SIC$long, "INV (L-H) within Chi quartile")
