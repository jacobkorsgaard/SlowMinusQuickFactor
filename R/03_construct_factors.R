## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026

## 3. Factor Construction #####################################################
## 3.0 Setup ##################################################################
cat("\014")
rm(list = ls())
graphics.off()

## Libraries
library(data.table)
library(lubridate)
library(zoo)

## Project paths
source(file.path("R", "00_config.R"))

## Load data (CRSP–Compustat monthly panel with adjustment cost estimates)
data <- readRDS("crsp_compustat_with_adj_cost.rds")
setDT(data)

## Basic sanity checks (optional but lightweight)
names(data)
uniqueN(data$permno)
uniqueN(data$gvkey)

## 3.1 Portfolio Formation Pipeline (June t → July t+1) #######################
## Timing convention:
## Monthly panel ffyear = t for July t through June t + 1.
## June of calendar year t is used to form portfolios/signals for holding year t.
## Since June itself belongs to ffyear = t - 1 in the monthly panel convention,
## June formation records are mapped using ffyear + 1L when constructing
## holding-period signals and portfolio assignments.
## 3.1.1 Timing, Sample Definition, and Eligibility ###########################
## Extract calendar components
data[, year := year(mdate)]
data[, month := month(mdate)]

## Fama–French portfolio year (ffyear = t for holding July t → June t+1)
data[, ffyear := fifelse(month >= 7, year, year - 1)]

## FF eligibility flags at formation (June t only)
june_eligible <-
  data[
    month == 6 &
      comp_ok == TRUE &
      has_dec == TRUE &
      has_june == TRUE
  ]

## 3.1.2 FF-style Annual Signals ##############################################
## Operating profitability (annual; fiscal year t-1 applied at June t)
op_ff <-
  data[
    !is.na(gvkey) & !is.na(fyear),
    .(
      revt     = revt[which.max(datadate)],
      cogs     = cogs[which.max(datadate)],
      xsga     = xsga[which.max(datadate)],
      xint     = xint[which.max(datadate)],
      be_clean = be_clean[which.max(datadate)]
    ),
    by = .(gvkey, fyear)
  ]


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

## Investment (asset growth; fiscal year t-1 applied at June t)
inv_ff <-
  data[
    !is.na(gvkey) & !is.na(fyear),
    .(
      at = at[which.max(datadate)]
    ),
    by = .(gvkey, fyear)
  ]

setorder(inv_ff, gvkey, fyear)
inv_ff[, at_lag := shift(at), by = gvkey]
inv_ff[, inv_ff :=
  fifelse(
    !is.na(at_lag) & at_lag > 0,
    (at - at_lag) / at_lag,
    NA_real_
  )]

## Align Compustat fiscal year to June formation year
inv_ff[, fyear := fyear + 1L]

## 3.1.3 June-t Formation Table ###############################################
## June t signals are attached to ffyear = t (holding July t → June t+1)
june_form <-
  data[
    month == 6 &
      has_dec == TRUE &
      has_june == TRUE,
    .(
      permno,
      gvkey,
      ffyear = ffyear + 1L,
      exchcd,
      me_june = me_clean,
      ppe_june = ppegt_lag,
      bm = bm_ff
    )
  ]

## Merge operating profitability (use fiscal year t-1 => fyear == ffyear-1)
## Implemented by merging on (gvkey, fyear) where fyear is stored in june_form as ffyear
## Merge operating profitability (FF timing)
june_form <-
  merge(
    june_form,
    op_ff[, .(gvkey, fyear, op_ff)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Merge investment (FF timing)
june_form <-
  merge(
    june_form,
    inv_ff[, .(gvkey, fyear, inv_ff)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Rename to match FF notation
setnames(june_form,
  old = c("op_ff", "inv_ff"),
  new = c("op", "inv")
)

## Enforce one row per permno–ffyear after merges
setorder(june_form, permno, ffyear)
june_form <- unique(june_form, by = c("permno", "ffyear"))


## 3.1.4 NYSE Breakpoints #####################################################
signals_nyse <- june_form[exchcd == 1]

## 30/70 breakpoints using quantiles (stable; avoids index edge cases)
ff_30_70_breaks <- function(x) {
  qs <- quantile(x, probs = c(0.3, 0.7), type = 1)
  c(low = as.numeric(qs[1]), high = as.numeric(qs[2]))
}

size_bp <-
  signals_nyse[, .(
    size_median = median(me_june, na.rm = TRUE)
  ), by = ffyear]

bm_bp <-
  signals_nyse[
    !is.na(bm),
    {
      b <- ff_30_70_breaks(bm)
      .(bm30 = b["low"], bm70 = b["high"])
    },
    by = ffyear
  ]

op_bp <-
  signals_nyse[
    !is.na(op),
    {
      b <- ff_30_70_breaks(op)
      .(op30 = b["low"], op70 = b["high"])
    },
    by = ffyear
  ]

inv_bp <-
  signals_nyse[
    !is.na(inv),
    {
      b <- ff_30_70_breaks(inv)
      .(inv30 = b["low"], inv70 = b["high"])
    },
    by = ffyear
  ]

june_form <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
    list(june_form, size_bp, bm_bp, op_bp, inv_bp)
  )

## 3.1.5 Portfolio Assignment #################################################
june_form[, size := fifelse(me_june <= size_median, "S", "B")]

june_form[, value :=
  fcase(
    bm <= bm30, "L",
    bm > bm70, "H",
    default = "M"
  )]

june_form[, prof :=
  fcase(
    op <= op30, "W",
    op > op70, "R",
    default = "M"
  )]

june_form[, invest :=
  fcase(
    inv <= inv30, "C",
    inv > inv70, "A",
    default = "M"
  )]

## Merge formation labels and frozen June weights into the monthly panel
## This applies June t info to ALL months in ffyear = t (July t .. June t+1)
data <-
  merge(
    data,
    june_form[, .(permno, ffyear, size, value, prof, invest, me_june, ppe_june)],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

## Time varying weights
setorder(data, permno, mdate)
data[, me_ff_w := shift(me_clean), by = permno]
data[me_ff_w <= 0, me_ff_w := NA_real_]

## 3.1.6 Monthly Portfolio Returns ############################################
## Require weights and returns
data <- data[!is.na(me_ff_w) & me_ff_w > 0]
# data = data[!is.na(me_june_w) & me_june_w > 0]

## 3.2 FF5 Factor Construction ################################################
## MKT is constructed as value-weighted market return minus RF
## SMB/HML/RMW/CMA are long–short factors (dollar neutral portfolios)

## 3.2.1 Market Factor ########################################################
## Risk-free rate by month
rf_m <- unique(data[!is.na(RF), .(mdate, RF)])
setkey(rf_m, mdate)

## Value-weighted market return using fixed June weights
mkt_vw <-
  data[, .(
    MKT_vw = weighted.mean(retadj, me_ff_w, na.rm = TRUE)
  ), by = mdate]
setkey(mkt_vw, mdate)

## Market excess return
mkt_rf <-
  rf_m[mkt_vw][, .(
    mdate,
    MKT = MKT_vw - RF
  )]

## 3.2.2 HML (High Minus Low) ################################################
hml <-
  data[!is.na(size) & !is.na(value) & value %in% c("H", "L"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, value)
  ][, .(
    HML = mean(ret[value == "H"]) - mean(ret[value == "L"])
  ), by = mdate]

## 3.2.3 SMB (Small Minus Big) ##############################################
## FF5 SMB is the average of SMB_HML, SMB_RMW, SMB_CMA

## SMB from size–value portfolios
smb_hml <-
  data[!is.na(size) & !is.na(value) & value %in% c("L", "M", "H"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, value)
  ]
smb_hml <-
  smb_hml[, .(
    SMB_HML = mean(ret[size == "S"]) - mean(ret[size == "B"])
  ), by = mdate]

## SMB from size–profitability portfolios
smb_rmw <-
  data[!is.na(size) & !is.na(prof) & prof %in% c("W", "M", "R"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, prof)
  ]
smb_rmw <-
  smb_rmw[, .(
    SMB_RMW = mean(ret[size == "S"]) - mean(ret[size == "B"])
  ), by = mdate]

## SMB from size–investment portfolios
smb_cma <-
  data[!is.na(size) & !is.na(invest) & invest %in% c("C", "M", "A"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, invest)
  ]
smb_cma <-
  smb_cma[, .(
    SMB_CMA = mean(ret[size == "S"]) - mean(ret[size == "B"])
  ), by = mdate]

## Combine SMB legs and average
smb <-
  Reduce(
    function(x, y) merge(x, y, by = "mdate", all = TRUE),
    list(smb_hml, smb_rmw, smb_cma)
  )
smb[, SMB := rowMeans(.SD, na.rm = TRUE),
  .SDcols = c("SMB_HML", "SMB_RMW", "SMB_CMA")
]
smb <- smb[, .(mdate, SMB)]

## 3.2.4 RMW (Robust Minus Weak) #############################################
rmw <-
  data[!is.na(size) & !is.na(prof) & prof %in% c("R", "W"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, prof)
  ][, .(
    RMW =
      0.5 * (ret[size == "S" & prof == "R"] - ret[size == "S" & prof == "W"]) +
        0.5 * (ret[size == "B" & prof == "R"] - ret[size == "B" & prof == "W"])
  ), by = mdate]

## 3.2.5 CMA (Conservative Minus Aggressive) ##################################
cma <-
  data[!is.na(size) & !is.na(invest) & invest %in% c("C", "A"),
    .(ret = weighted.mean(retadj, me_ff_w, na.rm = TRUE)),
    by = .(mdate, size, invest)
  ][, .(
    CMA =
      0.5 * (ret[size == "S" & invest == "C"] - ret[size == "S" & invest == "A"]) +
        0.5 * (ret[size == "B" & invest == "C"] - ret[size == "B" & invest == "A"])
  ), by = mdate]

## 3.2.6 Save FF5 ############################################################
ff5 <-
  Reduce(
    function(x, y) merge(x, y, by = "mdate", all = TRUE),
    list(mkt_rf, smb, hml, rmw, cma)
  )

saveRDS(ff5, "ff5_factors_internal.rds")

### 3.3 SMQ Factor Construction ################################################
### Goal: replicate FF factor mechanics exactly (timing + dynamic weights)
###   (i) Sorts use June formation information
###   (ii) Characteristic used for sorting is t-1 (lagged) to avoid look-ahead
###   (iii) Portfolio returns use continuously updated (dynamic) weights (same as FF legs)

### 3.3.1 Add industry identifiers to june_form ###############################
june_ids <-
  data[
    month == 6 &
      has_dec == TRUE &
      has_june == TRUE,
    .(
      permno,
      ffyear = ffyear + 1L,
      sic    = sic,
      sic2   = sic2
    )
  ]
setorder(june_ids, permno, ffyear)
june_ids <- unique(june_ids, by = c("permno", "ffyear"))

june_form <-
  merge(
    june_form,
    june_ids,
    by = c("permno", "ffyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

### 3.3.2 Build chi signals (t-1) aligned to June formation ###################
## Firm-level chi: (gvkey, fyear) then lag -> align like OP/INV (fyear + 1)
## Firm-level chi: mirror OP/INV exactly
chi_firm_ff <-
  data[
    !is.na(gvkey) & !is.na(fyear),
    .(
      chi_firm = chi_firm[which.max(datadate)]
    ),
    by = .(gvkey, fyear)
  ]
setorder(chi_firm_ff, gvkey, fyear)

## Lag one fiscal year
chi_firm_ff[, chi_firm_lag := shift(chi_firm, 1L), by = gvkey]

## Align to June formation year (fyear t-1 → June t)
chi_firm_ff[, fyear := fyear + 1L]

## Merge into June table
june_form <-
  merge(
    june_form,
    chi_firm_ff[, .(gvkey, fyear, chi_firm_lag)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Subindustry-level chi: (gsubind, fyear) then lag -> align (fyear + 1)
chi_sic_ff <-
  data[
    !is.na(sic) & !is.na(fyear),
    .(
      chi_sic = chi_sic[which.max(datadate)]
    ),
    by = .(sic, fyear)
  ]

setorder(chi_sic_ff, sic, fyear)
chi_sic_ff[, fyear := fyear + 1L]

june_form <-
  merge(
    june_form,
    chi_sic_ff[, .(sic, fyear, chi_sic)],
    by.x = c("sic", "ffyear"),
    by.y = c("sic", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Industry-level chi: (gind, fyear) then lag -> align (fyear + 1)
chi_sic2_ff <-
  data[
    !is.na(sic2) & !is.na(fyear),
    .(
      chi_sic2 = chi_sic2[which.max(datadate)]
    ),
    by = .(sic2, fyear)
  ]

setorder(chi_sic2_ff, sic2, fyear)
chi_sic2_ff[, chi_sic2_lag := shift(chi_sic2, 1L), by = sic2]
chi_sic2_ff[, fyear := fyear + 1L]

june_form <-
  merge(
    june_form,
    chi_sic2_ff[, .(sic2, fyear, chi_sic2_lag)],
    by.x = c("sic2", "ffyear"),
    by.y = c("sic2", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Enforce one row per permno–ffyear after chi merges (mirrors your pattern)
setorder(june_form, permno, ffyear)
june_form <- unique(june_form, by = c("permno", "ffyear"))

### 3.3.3 NYSE breakpoints for chi (30/70) ####################################
signals_nyse <- june_form[exchcd == 1]

chi_firm_bp <-
  signals_nyse[
    !is.na(chi_firm_lag),
    {
      b <- ff_30_70_breaks(chi_firm_lag)
      .(chi_firm_30 = b["low"], chi_firm_70 = b["high"])
    },
    by = ffyear
  ]

chi_sic_bp <-
  signals_nyse[
    !is.na(chi_sic),
    {
      b <- ff_30_70_breaks(chi_sic)
      .(chi_sic_30 = b["low"], chi_sic_70 = b["high"])
    },
    by = ffyear
  ]

chi_sic2_bp <-
  signals_nyse[
    !is.na(chi_sic2_lag),
    {
      b <- ff_30_70_breaks(chi_sic2_lag)
      .(chi_sic2_30 = b["low"], chi_sic2_70 = b["high"])
    },
    by = ffyear
  ]

### 3.3.4 Merge chi breakpoints + assign chi groups (Low/High; drop middle) ###
june_form <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
    list(june_form, chi_firm_bp, chi_sic_bp, chi_sic2_bp)
  )

june_form[, chi_firm_grp :=
  fcase(
    !is.na(chi_firm_lag) & chi_firm_lag <= chi_firm_30, "Low",
    !is.na(chi_firm_lag) & chi_firm_lag > chi_firm_70, "High",
    default = NA_character_
  )]

june_form[, chi_sic_grp :=
  fcase(
    !is.na(chi_sic) & chi_sic <= chi_sic_30, "Low",
    !is.na(chi_sic) & chi_sic > chi_sic_70, "High",
    default = NA_character_
  )]

june_form[, chi_sic2_grp :=
  fcase(
    !is.na(chi_sic2_lag) & chi_sic2_lag <= chi_sic2_30, "Low",
    !is.na(chi_sic2_lag) & chi_sic2_lag > chi_sic2_70, "High",
    default = NA_character_
  )]

### 3.3.5 Merge chi labels into monthly panel (like size/value/prof/invest) ####
data <-
  merge(
    data,
    june_form[, .(
      permno, ffyear,
      chi_firm_grp, chi_sic_grp, chi_sic2_grp
    )],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

### 3.3.6 SMQ helper (exactly mirrors HML/RMW/CMA mechanics) ##################
## - Uses existing FF size (S/B) from your pipeline
## - Uses dynamic weights (weight_var) exactly like you do with me_ff_w
## - Uses retadj and builds monthly long-short:
##     0.5*(S_High + B_High) - 0.5*(S_Low + B_Low)

make_SMQ <- function(dt, grp_var, factor_name, weight_var) {
  stopifnot(all(c("mdate", "size", "retadj") %in% names(dt)))
  stopifnot(grp_var %in% names(dt))
  stopifnot(weight_var %in% names(dt))

  port <-
    dt[
      !is.na(size) &
        size %in% c("S", "B") &
        !is.na(get(grp_var)) &
        get(grp_var) %in% c("Low", "High") &
        !is.na(get(weight_var)) &
        get(weight_var) > 0,
      {
        w <- get(weight_var) # <<< FIX: materialize weight
        .(ret = weighted.mean(retadj, w, na.rm = TRUE))
      },
      by = .(mdate, size, grp = get(grp_var))
    ]

  smq <-
    port[
      ,
      .(
        factor =
          0.5 * (ret[size == "S" & grp == "High"] +
            ret[size == "B" & grp == "High"]) -
            0.5 * (ret[size == "S" & grp == "Low"] +
              ret[size == "B" & grp == "Low"])
      ),
      by = mdate
    ]

  setnames(smq, "factor", factor_name)
  smq
}

### 3.3.7 Construct SMQ (Value-weighted, dynamic ME weights) ##################
## Value Weighted (Continously Updated)
SMQ_firm_VW <- make_SMQ(data, "chi_firm_grp", "SMQ_firm_VW", "me_ff_w")
SMQ_sic_VW <- make_SMQ(data, "chi_sic_grp", "SMQ_sic_VW", "me_ff_w")
SMQ_sic2_VW <- make_SMQ(data, "chi_sic2_grp", "SMQ_sic2_VW", "me_ff_w")

## Capital Weighted (Yearly Fixed)
data[, k_ff_w := ppe_june]
data[k_ff_w <= 0, k_ff_w := NA_real_]
SMQ_firm_K <- make_SMQ(data, "chi_firm_grp", "SMQ_firm_K", "k_ff_w")
SMQ_sic_K <- make_SMQ(data, "chi_sic_grp", "SMQ_sic_K", "k_ff_w")
SMQ_sic2_K <- make_SMQ(data, "chi_sic2_grp", "SMQ_sic2_K", "k_ff_w")

### 3.3.9 Save SMQ panel ######################################################
# smq_list = list(SMQ_firm_VW, SMQ_sub_VW, SMQ_ind_VW,
#                 SMQ_firm_K,  SMQ_sub_K,  SMQ_ind_K)
# lapply(smq_list, setkey, mdate)
#
# SMQ = Reduce(function(x, y) merge(x, y, by = "mdate", all = TRUE), smq_list)
#
# ## Restrict to main sample window used in the thesis (keep your window)
# SMQ_main = SMQ[mdate >= as.IDate("1981-07-01") & mdate <= as.IDate("2024-12-01")]
#
# saveRDS(SMQ_main, "SMQ_factors.rds")

### 3.3.X Expanding Chi Signals (t-1) #########################################

## Firm-level expanding chi
chi_firm_exp_ff <-
  data[
    !is.na(gvkey) & !is.na(fyear),
    .(
      chi_firm_exp = chi_firm_exp[which.max(datadate)]
    ),
    by = .(gvkey, fyear)
  ]

setorder(chi_firm_exp_ff, gvkey, fyear)
chi_firm_exp_ff[, chi_firm_exp_lag := shift(chi_firm_exp, 1L), by = gvkey]
chi_firm_exp_ff[, fyear := fyear + 1L]

june_form <-
  merge(
    june_form,
    chi_firm_exp_ff[, .(gvkey, fyear, chi_firm_exp_lag)],
    by.x = c("gvkey", "ffyear"),
    by.y = c("gvkey", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Industry-level expanding chi
chi_sic_exp_ff <-
  data[
    !is.na(sic) & !is.na(fyear),
    .(
      chi_sic_exp = chi_sic_exp[which.max(datadate)]
    ),
    by = .(sic, fyear)
  ]

setorder(chi_sic_exp_ff, sic, fyear)
chi_sic_exp_ff[, chi_sic_exp_lag := shift(chi_sic_exp, 1L), by = sic]
chi_sic_exp_ff[, fyear := fyear + 1L]

june_form <-
  merge(
    june_form,
    chi_sic_exp_ff[, .(sic, fyear, chi_sic_exp_lag)],
    by.x = c("sic", "ffyear"),
    by.y = c("sic", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

chi_sic2_exp_ff <-
  data[
    !is.na(sic2) & !is.na(fyear),
    .(
      chi_sic2_exp = chi_sic2_exp[which.max(datadate)]
    ),
    by = .(sic2, fyear)
  ]

setorder(chi_sic2_exp_ff, sic2, fyear)
chi_sic2_exp_ff[, chi_sic2_exp_lag := shift(chi_sic2_exp, 1L), by = sic2]
chi_sic2_exp_ff[, fyear := fyear + 1L]

june_form <-
  merge(
    june_form,
    chi_sic2_exp_ff[, .(sic2, fyear, chi_sic2_exp_lag)],
    by.x = c("sic2", "ffyear"),
    by.y = c("sic2", "fyear"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )

## Enforce uniqueness (same as baseline SMQs)
setorder(june_form, permno, ffyear)
june_form <- unique(june_form, by = c("permno", "ffyear"))

### Expanding chi NYSE breakpoints ############################################

signals_nyse <- june_form[exchcd == 1]

chi_firm_exp_bp <-
  signals_nyse[
    !is.na(chi_firm_exp_lag),
    {
      b <- ff_30_70_breaks(chi_firm_exp_lag)
      .(chi_firm_exp_30 = b["low"], chi_firm_exp_70 = b["high"])
    },
    by = ffyear
  ]

chi_sic_exp_bp <-
  signals_nyse[
    !is.na(chi_sic_exp_lag),
    {
      b <- ff_30_70_breaks(chi_sic_exp_lag)
      .(chi_sic_exp_30 = b["low"], chi_sic_exp_70 = b["high"])
    },
    by = ffyear
  ]

chi_sic2_exp_bp <-
  signals_nyse[
    !is.na(chi_sic2_exp_lag),
    {
      b <- ff_30_70_breaks(chi_sic2_exp_lag)
      .(chi_sic2_exp_30 = b["low"], chi_sic2_exp_70 = b["high"])
    },
    by = ffyear
  ]

### Expanding chi group assignment ############################################

june_form <-
  Reduce(
    function(x, y) merge(x, y, by = "ffyear", all.x = TRUE),
    list(
      june_form,
      chi_firm_exp_bp,
      chi_sic_exp_bp,
      chi_sic2_exp_bp
    )
  )

june_form[, chi_firm_exp_grp :=
  fcase(
    chi_firm_exp_lag <= chi_firm_exp_30, "Low",
    chi_firm_exp_lag > chi_firm_exp_70, "High",
    default = NA_character_
  )]

june_form[, chi_sic_exp_grp :=
  fcase(
    chi_sic_exp_lag <= chi_sic_exp_30, "Low",
    chi_sic_exp_lag > chi_sic_exp_70, "High",
    default = NA_character_
  )]

june_form[, chi_sic2_exp_grp :=
  fcase(
    chi_sic2_exp_lag <= chi_sic2_exp_30, "Low",
    chi_sic2_exp_lag > chi_sic2_exp_70, "High",
    default = NA_character_
  )]

data <-
  merge(
    data,
    june_form[, .(
      permno, ffyear,
      chi_firm_exp_grp,
      chi_sic_exp_grp,
      chi_sic2_exp_grp
    )],
    by = c("permno", "ffyear"),
    all.x = TRUE
  )

### Expanding-chi SMQs ########################################################
## Value-weighted (dynamic ME)
SMQ_firm_exp_VW <- make_SMQ(data, "chi_firm_exp_grp", "SMQ_firm_exp_VW", "me_ff_w")
SMQ_sic_exp_VW <- make_SMQ(data, "chi_sic_exp_grp", "SMQ_sic_exp_VW", "me_ff_w")
SMQ_sic2_exp_VW <- make_SMQ(data, "chi_sic2_exp_grp", "SMQ_sic2_exp_VW", "me_ff_w")

## Capital-weighted (June PPE)
SMQ_firm_exp_K <- make_SMQ(data, "chi_firm_exp_grp", "SMQ_firm_exp_K", "k_ff_w")
SMQ_sic_exp_K <- make_SMQ(data, "chi_sic_exp_grp", "SMQ_sic_exp_K", "k_ff_w")
SMQ_sic2_exp_K <- make_SMQ(data, "chi_sic2_exp_grp", "SMQ_sic2_exp_K", "k_ff_w")

### Update final SMQ panel #####################################################
smq_list <-
  list(
    SMQ_firm_VW, SMQ_sic_VW, SMQ_sic2_VW,
    SMQ_firm_K, SMQ_sic_K, SMQ_sic2_K,
    SMQ_firm_exp_VW, SMQ_sic_exp_VW, SMQ_sic2_exp_VW,
    SMQ_firm_exp_K, SMQ_sic_exp_K, SMQ_sic2_exp_K
  )

lapply(smq_list, setkey, mdate)

SMQ <-
  Reduce(function(x, y) merge(x, y, by = "mdate", all = TRUE), smq_list)

SMQ_main <-
  SMQ[mdate >= as.IDate("1981-07-01") & mdate <= as.IDate("2024-12-01")]

saveRDS(SMQ_main, "SMQ_factors.rds")
