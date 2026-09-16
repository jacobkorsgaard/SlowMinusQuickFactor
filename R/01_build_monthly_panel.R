## Firm Frictions, Saddle Paths and Risk Premia in the Cross-Section
## Authors: Jacob Korsgaard and Axel Emil Ulvemann
## Supervisor: Niels Joachim Gormsen
## Master Thesis 2026

## 1. Data Preparation ########################################################
## 1.1 Setup ##################################################################
# Clear Everything
cat("\014")
rm(list = ls())
graphics.off()

## Libraries
library(data.table)
library(lubridate)

## Project paths
source(file.path("R", "00_config.R"))

## 1.2 Load Raw Data ###########################################################
crsp <- fread(file.path(RAW_DATA_DIR, "CRSP-1950_01-2024_12.csv"))
cstat <- fread(file.path(RAW_DATA_DIR, "Compustat-1950_06 - 2024_12.csv"))
ccm <- fread(file.path(RAW_DATA_DIR, "Merger CRISPxCompustat.csv"))
rf <- fread(file.path(RAW_DATA_DIR, "CRSP - Riskfree - 1950_01_01 - 2024_12_31.csv"))

## 1.3 Standardize Variable Names #############################################
# Check Column Names
names(crsp)
names(cstat)
names(ccm)

# Function to clean column names
clean_names <- function(nm) {
  nm <- as.character(nm)
  nm <- gsub('^\\s*"(.*)"\\s*$', "\\1", nm)
  nm <- gsub("\\\\[a-zA-Z0-9]+\\s*", "", nm)
  nm <- gsub("\\\\+$", "", nm)
  nm <- gsub("[^A-Za-z0-9_]+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  tolower(nm)
}

# Apply Function
crsp <- copy(crsp)
setnames(crsp, clean_names(names(crsp)))
cstat <- copy(cstat)
setnames(cstat, clean_names(names(cstat)))
ccm <- copy(ccm)
setnames(ccm, clean_names(names(ccm)))

## 1.4 Parse Date Fields #######################################################
# Function to Parse CCM Date Fields
parse_ccm_yyyymmdd <- function(x) {
  if (is.numeric(x)) x <- format(x, scientific = FALSE, trim = TRUE)
  x <- trimws(as.character(x))
  out <- rep(as.Date(NA), length(x))
  idx <- grepl("^[0-9]{8}$", x)
  out[idx] <- as.Date(x[idx], format = "%Y%m%d")
  idx <- is.na(out) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  out[idx] <- as.Date(x[idx], format = "%Y-%m-%d")
  out
}

## 1.5 CRSP Variable Types and Monthly Timing #################################
# Enforce CRSP variable types and construct monthly time identifiers
crsp[, permno := as.integer(permno)]
crsp[, date := as.IDate(date)]
crsp[, mdate := as.IDate(floor_date(as.Date(date), "month"))]
crsp[, yyyymm := year(as.Date(mdate)) * 100L + month(as.Date(mdate))]

# Coerce CRSP return and price fields to numeric
num_cols_crsp <- intersect(c("ret", "dlret", "prc"), names(crsp))
for (v in num_cols_crsp) crsp[, (v) := as.numeric(gsub(",", "", as.character(get(v))))]

# Clean and coerce shares outstanding to numeric
crsp[, shrout_raw := as.character(shrout)]
crsp[, shrout_clean := trimws(gsub("\\\\", "", shrout_raw))]
crsp[shrout_clean == "", shrout_clean := NA_character_]
crsp[, shrout := as.numeric(gsub("[^0-9\\.\\-]", "", shrout_clean))]

# Enforce CRSP exchange and share class codes as integers
crsp[, exchcd := as.integer(exchcd)]
crsp[, shrcd := as.integer(shrcd)]

## 1.6 Compustat Variable Types ################################################
# Compustat types
cstat[, gvkey := as.character(gvkey)]
cstat[, datadate := as.IDate(datadate)]
cstat[, fyear := as.integer(fyear)]

## 1.7 CCM Link Table Preparation #############################################
# CCM types and dates
ccm[, gvkey := as.character(gvkey)]
ccm[, lpermno := as.integer(lpermno)]
ccm[, linkdt := parse_ccm_yyyymmdd(linkdt)]
ccm[, linkenddt := parse_ccm_yyyymmdd(linkenddt)]
ccm[is.na(linkenddt), linkenddt := as.Date("9999-12-31")]

# Keep standard CRSP–Compustat links
ccm <- ccm[linktype %in% c("LU", "LC", "LN", "LS")]
ccm <- ccm[linkprim %in% c("P", "C")]

## 1.8 Initialize CRSP Working Dataset ########################################
# Working dataset for CRSP–Compustat integration
data <- copy(crsp)

## 1.9 CRSP–Compustat Linking ##################################################
# Prepare CCM for linking
ccm <- ccm[, .(gvkey, permno = lpermno, linkdt, linkenddt, linktype, linkprim)]
setorder(data, permno, mdate)
setorder(ccm, permno, linkdt, linkenddt)

setkey(data, permno)
setkey(ccm, permno, linkdt, linkenddt)

# Monthly CRSP–Compustat link
data <- ccm[
  data,
  on = .(permno, linkdt <= mdate, linkenddt >= mdate),
  nomatch = 0L,
  allow.cartesian = TRUE,
  .(permno,
    mdate  = i.mdate,
    yyyymm = i.yyyymm,
    exchcd = i.exchcd,
    shrcd  = i.shrcd,
    prc    = i.prc,
    shrout = i.shrout,
    ret    = i.ret,
    dlret  = i.dlret,
    gvkey,
    linkprim,
    linktype,
    linkdt
  )
]

# Rank links
data[, linkprim_rank := fifelse(
  linkprim == "P", 1L,
  fifelse(linkprim == "C", 2L, 9L)
)]
data[, linktype_rank := fcase(
  linktype == "LU", 1L,
  linktype == "LC", 2L,
  linktype == "LN", 3L,
  linktype == "LS", 4L,
  default = 9L
)]

# Prefer most recent valid link
data[, linkdt_num := as.integer(as.Date(linkdt))]
setorder(data, permno, mdate, linkprim_rank, linktype_rank, -linkdt_num)

# Keep one gvkey per permno–month
data <- data[, .SD[1], by = .(permno, mdate, yyyymm)]
setkey(data, permno, mdate)

## 1.10 Compustat Numeric Variable Cleaning ###################################
# Compustat numeric variables to retain/clean
comp_num <- intersect(
  c(
    "at", "act", "ppegt", "dp", "capx", "pstk", "pstkl", "pstkrv", "seq",
    "txditc", "cogs", "revt", "xint", "xsga", "dlc", "dltt", "prcc_f", "csho"
  ),
  names(cstat)
)

# Helper to coerce Compustat variables to numeric
to_num <- function(x) {
  x0 <- as.character(x)
  x0 <- gsub("\\\\", "", x0)
  x0 <- gsub(",", "", x0)
  x0 <- gsub("[^0-9\\.\\-]", "", x0)
  x0[x0 %in% c("", ".", "-", "-.")] <- NA
  as.numeric(x0)
}

# Apply to Compustat numeric variables
for (v in comp_num) cstat[, (v) := to_num(get(v))]

## 1.11 July Convention Timing #################################################
# Compustat industry classification variables
comp_ind <- intersect(c("gsector", "gind", "gsubind", "sic"), names(cstat))

# July convention availability window
cstat[, start_mdate := as.IDate(make_date(year(as.Date(datadate)) + 1L, 7L, 1L))]
cstat[, end_mdate := as.IDate(make_date(year(as.Date(datadate)) + 2L, 6L, 1L))]

# Keep required Compustat fields
comp_keep <- intersect(
  c("gvkey", "datadate", "fyear", "start_mdate", "end_mdate", comp_num, comp_ind),
  names(cstat)
)
cstat <- cstat[, ..comp_keep]

## 1.12 Merge Compustat into CRSP ##############################################
# Prepare keys for lagged merge
data[, mdate_i := mdate]
setkey(cstat, gvkey, start_mdate, end_mdate)
setkey(data, gvkey)

# Merge Compustat into CRSP using July convention
data <- cstat[
  data,
  on = .(gvkey, start_mdate <= mdate_i, end_mdate >= mdate_i),
  allow.cartesian = TRUE
]

# Tie-break: keep most recent accounting data
data[, datadate_num := as.integer(datadate)]
setorder(data, permno, mdate, -datadate_num)
data <- data[, .SD[1], by = .(permno, mdate)]

# Final cleanup (link ranks and helper keys)
data[, `:=`(
  datadate_num = NULL,
  mdate_i = NULL,
  linkprim_rank = NULL,
  linktype_rank = NULL,
  linkdt_num = NULL
)]

## 1.13 Define FF Universe ####################################################
# FF universe: common stocks on NYSE, AMEX, NASDAQ
data <- data[exchcd %in% c(1L, 2L, 3L)]
data <- data[shrcd %in% c(10L, 11L)]

## 1.14 Merge Risk-Free Rate ##################################################
setnames(rf, tolower(gsub("[^A-Za-z0-9_]", "_", names(rf))))
rf[, mcaldt := as.IDate(mcaldt)]
rf[, mdate := as.IDate(floor_date(as.Date(mcaldt), "month"))]

# Convert annualized yield (%) to monthly simple return (assume 30-day month)
rf[, rf := exp((tmytm / 100) * (30 / 365)) - 1]
rf <- rf[, .(RF = mean(rf, na.rm = TRUE)), by = mdate]

setkey(rf, mdate)
setkey(data, mdate)
data <- rf[data]

## 1.15 Final Variable Retention ##############################################
## Create sic2
data[, sic := as.integer(sic)]
data[, sic2 := fifelse(!is.na(sic), sic %/% 100L, NA_integer_)]

vars_keep <- c(
  # Identifiers & time
  "permno", "gvkey", "mdate", "yyyymm", "fyear", "datadate",

  # Industry classification
  "gsector", "gind", "gsubind", "sic", "sic2",

  # Returns & prices (raw only; economics in Section 2)
  "ret", "dlret", "RF", "prc", "shrout",

  # Raw Compustat inputs
  "at", "act", "ppegt", "capx",
  "dlc", "dltt", "prcc_f", "seq", "txditc", "csho",
  "pstkrv", "pstkl", "pstk",
  "revt", "cogs", "xsga", "xint",

  # Controls
  "exchcd", "shrcd"
)

data <- data[, intersect(vars_keep, names(data)), with = FALSE]

## 1.16 Save CRSP–Compustat Monthly Panel ######################################
saveRDS(data, "crsp_compustat_monthly_panel_1950_2024.rds")
