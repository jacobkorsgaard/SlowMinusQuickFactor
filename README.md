# Firm Frictions, Adjustment Dynamics, and Cross-Sectional Returns

This repository contains the R code accompanying the 2026 Copenhagen Business
School master's thesis *Firm Frictions, Adjustment Dynamics, and Cross-Sectional
Returns* by Jacob Korsgaard and Axel Ulveman.

The study examines whether the Fama–French five-factor model captures economic
mechanisms associated with the q-theory of investment. It estimates firms' speed
of adjustment from rolling investment–q regressions and constructs a
slow-minus-quick (SMQ) factor from the resulting adjustment-cost estimates.

## Main result

SMQ does not earn a clear unconditional risk premium, but it earns a positive
and statistically significant alpha relative to the Fama–French five-factor
model. The pricing evidence is concentrated primarily among small firms. Adding
SMQ to the five-factor model produces only limited improvements in overall
asset-pricing performance.

The complete thesis is available in
[`final-thesis/`](final-thesis/Firm%20Frictions,%20Adjustment%20Dynamics,%20and%20Cross-Sectional%20Returns.pdf).

## Repository structure

```text
.
├── R/                  # Ordered data-construction and analysis scripts
├── data/
│   ├── raw/            # User-supplied source data (not tracked by Git)
│   └── interim/        # Generated R objects (not tracked by Git)
├── final-thesis/       # Final thesis PDF
├── results/
│   ├── figures/        # Generated figures (not tracked by Git)
│   └── tables/         # Reserved for exported tables (not tracked by Git)
└── run_all.R           # Sequential pipeline runner
```

## Code workflow

| Script | Purpose |
|---|---|
| `00_config.R` | Defines portable project paths and output directories |
| `01_build_monthly_panel.R` | Cleans and merges CRSP, Compustat, CCM, and risk-free data |
| `02_estimate_adjustment_costs.R` | Constructs q-theory variables and rolling adjustment-cost estimates |
| `03_construct_factors.R` | Replicates FF5 and constructs SMQ variants |
| `04_construct_test_portfolios.R` | Forms one-, two-, and three-dimensional test portfolios |
| `05_build_analysis_panel.R` | Merges factors and test assets into the analysis panel |
| `06_summary_statistics.R` | Produces characteristic and adjustment-cost diagnostics |
| `07_analyze_smq_strategy.R` | Analyzes SMQ, adjustment-cost deciles, and mean–variance frontiers |
| `08_asset_pricing_tests.R` | Runs factor spanning, time-series, and GRS tests |
| `09_subperiod_analysis.R` | Studies subperiod evidence in the investment–q relationship |
| `10_robustness_checks.R` | Evaluates alternative adjustment-cost specifications and robustness |
| `11_appendix_validation.R` | Validates replicated factors and portfolios against French data |

Scripts are numbered in their intended execution order. The original empirical
definitions, sample filters, estimators, and portfolio-construction rules are
retained.

## Data

The data are deliberately excluded from this repository. CRSP and Compustat are
licensed products and cannot be redistributed here. Researchers must obtain
access independently and place the following files in `data/raw/`:

```text
CRSP-1950_01-2024_12.csv
Compustat-1950_06 - 2024_12.csv
Merger CRISPxCompustat.csv
CRSP - Riskfree - 1950_01_01 - 2024_12_31.csv
FF5 Download - 1963_07-2024_12.csv
kf_portfolios/25_Portfolios_5x5.csv
kf_portfolios/25_Portfolios_ME_OP_5x5.csv
kf_portfolios/25_Portfolios_ME_INV_5x5.csv
kf_portfolios/32_Portfolios_ME_OP_INV_2x4x4.csv
kf_portfolios/32_Portfolios_ME_BEME_INV_2x4x4.csv
kf_portfolios/32_Portfolios_ME_BEME_OP_2x4x4.csv
```

The Kenneth French comparison files can be downloaded from the Kenneth R.
French Data Library. Their local filenames must match the list above.

No raw data, generated `.rds` objects, tables, or figures are tracked by Git;
see [`.gitignore`](.gitignore).

## Requirements

The analysis requires R and these packages:

```r
install.packages(c(
  "data.table", "ggplot2", "lmtest", "lubridate", "moments",
  "sandwich", "xtable", "zoo"
))
```

The thesis computations were performed using monthly CRSP data through December
2024 and annual Compustat data covering the associated sample period.

## Reproducing the analysis

Run commands from the repository root. To execute the complete pipeline:

```sh
Rscript run_all.R
```

Individual stages can also be run separately, provided their upstream `.rds`
inputs already exist. For example:

```sh
Rscript R/03_construct_factors.R
```

Generated intermediate datasets and R table objects are written to
`data/interim/`; figures are written to `results/figures/`. The full pipeline is
computationally demanding, particularly the rolling adjustment-cost estimation
and robustness analysis.

## Verification

The complete pipeline was successfully executed after restructuring. Results
are consistent with the thesis conclusions, although individual estimates may
differ slightly because of data vintages, package versions, and subsequent code
cleanup.

## Citation

Korsgaard, Jacob, and Axel Ulveman (2026). *Firm Frictions, Adjustment Dynamics,
and Cross-Sectional Returns*. Master's thesis, Copenhagen Business School.
