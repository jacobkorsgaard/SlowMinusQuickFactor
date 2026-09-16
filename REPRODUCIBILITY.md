# Reproducibility status

## Verification performed

The complete pipeline was executed from the raw files currently stored in
`data/raw/`. All eleven active scripts completed successfully after the
repository cleanup.

The run produced:

- a 3,280,892-row CRSP–Compustat monthly panel covering January 1950 through
  December 2024;
- rolling firm-, four-digit SIC-, and two-digit SIC-level adjustment-cost
  estimates;
- internally constructed Fama–French factors and twelve SMQ variants;
- one-, two-, and three-dimensional test portfolios;
- summary statistics, factor-spanning regressions, GRS tests, robustness
  results, figures, and appendix validation outputs.

The replicated Fama–French factors closely track the Kenneth French benchmark.
For July 1982–December 2024, own-versus-benchmark correlations are approximately
1.000 for MKT, 0.995 for SMB, 0.962 for HML, 0.968 for RMW, and 0.956 for CMA.

## Difference from the submitted thesis

The clean run does not reproduce every printed number in the final thesis PDF.
Most notably, the current pipeline produces an FF5 alpha for the baseline SMQ
factor of 0.37% per month (Newey–West t-statistic 3.35), whereas Table 1 of the
submitted thesis reports 0.32% (t-statistic 3.19). The current SMQ series has 504
non-missing observations from July 1982 through June 2024, while the thesis table
is labelled July 1982–December 2024.

This difference was present when running the supplied code and data and was not
introduced by filename, path, formatting, or robustness-script repairs. It
indicates version drift between the archived thesis results and the final code or
input extracts. The repository therefore distinguishes between:

1. **Computational reproducibility:** the published pipeline runs end to end;
2. **Numerical reproduction:** exact agreement with every PDF table remains
   pending identification of the precise final-thesis factor/input snapshot.

No empirical definitions have been altered merely to force agreement with a
reported headline result.

## Repairs required for a clean run

The publication cleanup made the following non-substantive or correctness
repairs:

- replaced machine-specific working directories with project-relative paths;
- corrected packaged input filenames and Kenneth French benchmark paths;
- parsed the supplied RTF-wrapped Fama–French CSV correctly;
- aligned orthogonalized robustness-factor residuals with their original dates
  after complete-case regression filtering;
- supplied missing family-order, prefix-mapping, quartile-label, and display
  helpers in the final robustness-table blocks.

The row-alignment repair changes only previously failing robustness code. The
other repairs concern paths, labels, or output assembly.
