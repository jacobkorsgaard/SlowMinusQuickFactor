### 6. Summary Statistics and Chi Diagnostics ################################
cat("\014")
rm(list = ls())
graphics.off()

### 6.1 Setup #################################################################
library(data.table)
library(lubridate)
library(moments)

source(file.path("R", "00_config.R"))

START_DATE <- as.Date("1982-07-01")
END_DATE <- as.Date("2024-12-01")

### 6.2 Load Data #############################################################
panel <- readRDS("crsp_compustat_with_adj_cost.rds")
setDT(panel)
panel[, mdate := as.Date(mdate)]

### 6.3 Calendar Fields #######################################################
panel[, year := year(mdate)]
panel[, month := month(mdate)]
panel[, ffyear := fifelse(month >= 7, year, year - 1)]
panel[, sic := as.integer(sic)]

if (!"sic2" %in% names(panel)) {
  panel[, sic2 := fifelse(!is.na(sic), sic %/% 10L, NA_integer_)]
}

### 6.4 June Formation Sample #################################################
form_sample <-
  panel[
    month == 6 &
      has_dec == TRUE &
      has_june == TRUE &
      mdate >= START_DATE &
      mdate <= END_DATE
  ]

### 6.5 Helpers ###############################################################
winsor_minmax <- function(x, p = 0.01) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(c(wmin = NA_real_, wmax = NA_real_))
  }
  qs <- quantile(x, probs = c(p, 1 - p), na.rm = TRUE, type = 7)
  c(wmin = as.numeric(qs[1]), wmax = as.numeric(qs[2]))
}

winsor_1pct <- function(x) {
  x_num <- as.numeric(x)
  q <- quantile(x_num, probs = c(0.01, 0.99), na.rm = TRUE, type = 7)
  pmin(pmax(x_num, q[1]), q[2])
}

safe_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) {
    return(NA_real_)
  }
  skewness(x)
}

safe_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4) {
    return(NA_real_)
  }
  kurtosis(x)
}

make_summary_row <- function(x, var_name, add_shape = FALSE) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    out <- data.table(
      Variable = var_name,
      Mean = NA_real_,
      Median = NA_real_,
      SD = NA_real_,
      P10 = NA_real_,
      P25 = NA_real_,
      P75 = NA_real_,
      P90 = NA_real_,
      N = 0L
    )

    if (add_shape) {
      out[, `:=`(
        Skewness = NA_real_,
        Kurtosis = NA_real_
      )]
    }

    return(out)
  }

  out <- data.table(
    Variable = var_name,
    Mean = mean(x),
    Median = median(x),
    SD = sd(x),
    P10 = quantile(x, 0.10, na.rm = TRUE, type = 7),
    P25 = quantile(x, 0.25, na.rm = TRUE, type = 7),
    P75 = quantile(x, 0.75, na.rm = TRUE, type = 7),
    P90 = quantile(x, 0.90, na.rm = TRUE, type = 7),
    N = length(x)
  )

  if (add_shape) {
    out[, `:=`(
      Skewness = safe_skewness(x),
      Kurtosis = safe_kurtosis(x)
    )]
  }

  out
}

safe_min_year <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  min(x)
}

safe_max_year <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

### 6.6 Variables #############################################################
char_vars <- c(
  "bm_ff", "op_clean", "inv_annual", "ik_annual", "tobins_q",
  "chi_firm", "chi_sic", "chi_sic2"
)

char_labels <- c(
  bm_ff      = "BM",
  op_clean   = "OP",
  inv_annual = "INV",
  ik_annual  = "I/K",
  tobins_q   = "q",
  chi_firm   = "chi_firm",
  chi_sic    = "chi_sic",
  chi_sic2   = "chi_sic2"
)

chi_vars <- c("chi_firm", "chi_sic", "chi_sic2")

plot_vars <- c(
  "bm_ff", "op_clean", "inv_annual",
  "chi_firm", "chi_sic", "chi_sic2"
)

### 6.7 Characteristic Summary ################################################
characteristic_summary <-
  rbindlist(
    lapply(char_vars, function(v) {
      make_summary_row(form_sample[[v]], char_labels[[v]], add_shape = FALSE)
    }),
    fill = TRUE
  )

characteristic_summary_appendix <-
  rbindlist(
    lapply(char_vars, function(v) {
      make_summary_row(form_sample[[v]], char_labels[[v]], add_shape = TRUE)
    }),
    fill = TRUE
  )

cat("\n============================================================\n")
cat("CHARACTERISTIC SUMMARY\n")
cat("============================================================\n")
print(characteristic_summary)

cat("\n============================================================\n")
cat("CHARACTERISTIC SUMMARY APPENDIX\n")
cat("============================================================\n")
print(characteristic_summary_appendix)

### 6.8 Winsorized Characteristic Summary #####################################
for (v in char_vars) {
  form_sample[, paste0(v, "_w1") := winsor_1pct(get(v))]
}

characteristic_summary_winsorized <-
  rbindlist(
    lapply(char_vars, function(v) {
      make_summary_row(
        form_sample[[paste0(v, "_w1")]],
        char_labels[[v]],
        add_shape = TRUE
      )
    }),
    fill = TRUE
  )

cat("\n============================================================\n")
cat("WINSORIZED CHARACTERISTIC SUMMARY (1%)\n")
cat("============================================================\n")
print(characteristic_summary_winsorized)

### 6.9 Formation Sample Missingness ##########################################
formation_missingness <-
  rbindlist(
    lapply(char_vars, function(v) {
      x <- form_sample[[v]]
      data.table(
        Variable = char_labels[[v]],
        N_total = length(x),
        N_valid = sum(!is.na(x)),
        N_missing = sum(is.na(x)),
        Share_missing = mean(is.na(x))
      )
    }),
    fill = TRUE
  )

cat("\n============================================================\n")
cat("FORMATION SAMPLE MISSINGNESS\n")
cat("============================================================\n")
print(formation_missingness)

### 6.10 Chi Distribution Diagnostics #########################################
chi_distribution_table <-
  rbindlist(
    lapply(chi_vars, function(v) {
      x <- form_sample[[v]]
      x_clean <- x[is.finite(x)]
      wm <- winsor_minmax(x_clean, p = 0.01)

      data.table(
        Variable = v,
        Mean = mean(x_clean, na.rm = TRUE),
        Median = median(x_clean, na.rm = TRUE),
        SD = sd(x_clean, na.rm = TRUE),
        Skewness = safe_skewness(x_clean),
        Kurtosis = safe_kurtosis(x_clean),
        P10 = quantile(x_clean, 0.10, na.rm = TRUE, type = 7),
        P90 = quantile(x_clean, 0.90, na.rm = TRUE, type = 7),
        Winsor_P1 = wm["wmin"],
        Winsor_P99 = wm["wmax"],
        N = length(x_clean)
      )
    }),
    fill = TRUE
  )

cat("\n============================================================\n")
cat("CHI DISTRIBUTION DIAGNOSTICS\n")
cat("============================================================\n")
print(chi_distribution_table)

### 6.11 Chi Estimation Quality Diagnostics ###################################
chi_quality_table <-
  rbindlist(list(
    data.table(
      Level = "Firm",
      ID_Variable = "gvkey",
      Year_Variable = "fyear",
      N_obs = panel[!is.na(chi_firm), .N],
      N_unique_ids = uniqueN(panel[!is.na(chi_firm), gvkey]),
      N_unique_id_year = uniqueN(panel[!is.na(chi_firm), paste(gvkey, fyear)]),
      First_year = safe_min_year(panel[!is.na(chi_firm), fyear]),
      Last_year = safe_max_year(panel[!is.na(chi_firm), fyear])
    ),
    data.table(
      Level = "SIC",
      ID_Variable = "sic",
      Year_Variable = "fyear",
      N_obs = panel[!is.na(chi_sic), .N],
      N_unique_ids = uniqueN(panel[!is.na(chi_sic), sic]),
      N_unique_id_year = uniqueN(panel[!is.na(chi_sic), paste(sic, fyear)]),
      First_year = safe_min_year(panel[!is.na(chi_sic), fyear]),
      Last_year = safe_max_year(panel[!is.na(chi_sic), fyear])
    ),
    data.table(
      Level = "SIC2",
      ID_Variable = "sic2",
      Year_Variable = "fyear",
      N_obs = panel[!is.na(chi_sic2), .N],
      N_unique_ids = uniqueN(panel[!is.na(chi_sic2), sic2]),
      N_unique_id_year = uniqueN(panel[!is.na(chi_sic2), paste(sic2, fyear)]),
      First_year = safe_min_year(panel[!is.na(chi_sic2), fyear]),
      Last_year = safe_max_year(panel[!is.na(chi_sic2), fyear])
    )
  ), fill = TRUE)

cat("\n============================================================\n")
cat("CHI ESTIMATION QUALITY DIAGNOSTICS\n")
cat("============================================================\n")
print(chi_quality_table)

### 6.12 Density Diagnostics ##################################################
cat("\n============================================================\n")
cat("DENSITY DIAGNOSTICS: RAW VARIABLES\n")
cat("============================================================\n")

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (v in plot_vars) {
  x <- form_sample[[v]]
  x <- x[is.finite(x)]

  if (length(x) > 1) {
    plot(
      density(x, na.rm = TRUE),
      main = paste("Density:", v),
      xlab = v
    )
  } else {
    plot.new()
    title(main = paste("Density:", v))
    text(0.5, 0.5, "Insufficient data")
  }
}

par(mfrow = c(1, 1))

cat("\n============================================================\n")
cat("DENSITY DIAGNOSTICS: WINSORIZED AT 1%\n")
cat("============================================================\n")

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (v in plot_vars) {
  vw <- paste0(v, "_w1")
  x <- form_sample[[vw]]
  x <- x[is.finite(x)]

  if (length(x) > 1) {
    plot(
      density(x, na.rm = TRUE),
      main = paste("Density:", v, "(winsorized 1%)"),
      xlab = v
    )
  } else {
    plot.new()
    title(main = paste("Density:", v, "(winsorized 1%)"))
    text(0.5, 0.5, "Insufficient data")
  }
}

par(mfrow = c(1, 1))

### 6.13 Save Outputs #########################################################
saveRDS(characteristic_summary, "table_characteristic_summary.rds")
saveRDS(characteristic_summary_appendix, "table_characteristic_summary_appendix.rds")
saveRDS(characteristic_summary_winsorized, "table_characteristic_summary_winsorized_1pct.rds")
saveRDS(formation_missingness, "table_formation_sample_missingness.rds")
saveRDS(chi_distribution_table, "table_chi_distribution_diagnostics.rds")
saveRDS(chi_quality_table, "table_chi_quality_base.rds")
