#!/usr/bin/env Rscript

# Run each stage in a clean R session so that the original scripts' workspace
# cleanup cannot remove state needed by the runner.
scripts <- sprintf(
  "R/%02d_%s.R",
  1:11,
  c(
    "build_monthly_panel",
    "estimate_adjustment_costs",
    "construct_factors",
    "construct_test_portfolios",
    "build_analysis_panel",
    "summary_statistics",
    "analyze_smq_strategy",
    "asset_pricing_tests",
    "subperiod_analysis",
    "robustness_checks",
    "appendix_validation"
  )
)

for (script in scripts) {
  message("\nRunning ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) {
    stop("Pipeline failed in ", script, call. = FALSE)
  }
}

message("\nPipeline completed successfully.")
