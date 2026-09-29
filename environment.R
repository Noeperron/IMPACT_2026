# ============================================================================
# Purpose:      One-shot R package installer for the figure scripts. Installs the CRAN packages
#               used across Figure1-5 and SupplementalFigures, pinned to the versions the
#               published figures were generated with so a fresh install reproduces them
#               rather than tracking upstream releases.
# Inputs:       None (idempotent: a package already present at the pinned version is skipped).
# Outputs:      Installed R packages, listed in PACKAGES below.
# Dependencies: R 4.5+ with internet access to CRAN. Pinned installs additionally need the
#               'remotes' package, which is installed automatically if missing.
# ============================================================================

# Pinned to the environment used for the published figures (R 4.5.1).
PACKAGES <- c(
  tidyverse = "2.0.0",
  ggplot2   = "4.0.1",
  ggpubr    = "0.6.2",
  rstatix   = "0.7.3",
  ggrepel   = "0.9.6",
  plotrix   = "3.8.13",
  lme4      = "1.1.38",
  lmerTest  = "3.2.0"
)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

for (pkg in names(PACKAGES)) {
  want <- PACKAGES[[pkg]]
  have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
  if (identical(have, want)) {
    cat(sprintf("%-10s %s already installed\n", pkg, want))
  } else {
    cat(sprintf("%-10s installing %s (found: %s)\n", pkg, want,
                ifelse(is.na(have), "none", have)))
    remotes::install_version(pkg, version = want, repos = "https://cloud.r-project.org",
                             upgrade = "never")
  }
}

cat("\nAll packages installed at the pinned versions.\n")
