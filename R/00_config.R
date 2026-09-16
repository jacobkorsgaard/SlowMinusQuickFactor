# Shared project paths ---------------------------------------------------------
# Run scripts from the repository root. All generated data are written to
# data/interim, which is excluded from version control.

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!file.exists(file.path(project_root, "R", "00_config.R"))) {
  stop(
    "Run this script from the repository root (the directory containing R/).",
    call. = FALSE
  )
}

RAW_DATA_DIR <- file.path(project_root, "data", "raw")
INTERIM_DATA_DIR <- file.path(project_root, "data", "interim")
TABLES_DIR <- file.path(project_root, "results", "tables")
FIGURES_DIR <- file.path(project_root, "results", "figures")

dir.create(INTERIM_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# The original thesis scripts read and wrote intermediate objects in one flat
# working directory. Keeping data/interim as the working directory preserves
# that behavior without changing any empirical calculations.
setwd(INTERIM_DATA_DIR)
