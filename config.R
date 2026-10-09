# ============================================================================
# Purpose:      Shared configuration for R figure scripts: data, repo and figure-output paths, color palette, and the age+sex-adjusted Jonckheere-Terpstra helper.
# Inputs:       Environment variable expectations; the user edits SCRNA_DIR to the directory containing scRNAseq_IMPACT_Zenodo.h5ad.
# Outputs:      DATA_DIR, FIGURES_DIR, TABLES_DIR, COLORS, FONT, save_figure(), jt_test_residuals_age_sex() (sourced by every R figure script).
# Dependencies: R base; tryCatch; rstudioapi (optional, command-line Rscript falls back to getwd()).
# ============================================================================
# IMPACT configuration
# Edit these paths to match your local environment.

# Base directory containing scRNA-seq h5ad files (from GEO/Zenodo)
SCRNA_DIR <- "/path/to/scrnaseq_data"

# Base directory for this repository (auto-detected; works in RStudio AND command-line Rscript)
REPO_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) ""
)
if (is.null(REPO_DIR) || REPO_DIR == "") {
  REPO_DIR <- getwd()
}

# Data files (download from Zenodo: https://doi.org/10.5281/zenodo.18989222)
DATA_DIR <- file.path(REPO_DIR, "..", "data")

# Output directory for figures
FIGURES_DIR <- file.path(REPO_DIR, "..", "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# Output directory for generated tables and figure source data
TABLES_DIR <- file.path(REPO_DIR, "..", "tables")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# TYPOGRAPHY. Every figure in this paper is set in Arial. Arial is not
# redistributable and is usually absent on Linux, so if it is missing we register
# Liberation Sans under that name. The two are metrically identical (identical
# glyph advance widths), so text occupies exactly the same space and nothing
# reflows. Registering rather than relying on silent substitution matters for
# vector output: the SVG then records the family we asked for.
FONT <- "Arial"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  if (!(FONT %in% systemfonts::system_fonts()$family)) {
    .sub <- systemfonts::system_fonts()
    .sub <- .sub[.sub$family == "Liberation Sans", ]
    if (nrow(.sub) > 0) {
      .pick <- function(st) {
        hit <- .sub$path[.sub$style == st]
        if (length(hit)) hit[1] else .sub$path[.sub$style == "Regular"][1]
      }
      systemfonts::register_font(FONT, plain = .pick("Regular"), bold = .pick("Bold"),
                                 italic = .pick("Italic"), bolditalic = .pick("Bold Italic"))
    }
  }
}

# Write a ggplot to png, pdf and svg with Arial preserved.
#   png via ragg and pdf via cairo, because the base devices cannot resolve a
#   named font family and error out instead of substituting.
#   svg via svglite, which keeps text as text (editable), unlike the cairo SVG
#   device which converts glyphs to outlines. svglite records the RESOLVED family
#   name, so the file is post-processed to say Arial.
save_figure <- function(plot, basename, width, height, dpi = 300, svg = TRUE) {
  ggplot2::ggsave(file.path(FIGURES_DIR, paste0(basename, ".png")), plot,
                  width = width, height = height, units = "in", dpi = dpi,
                  device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL)
  ggplot2::ggsave(file.path(FIGURES_DIR, paste0(basename, ".pdf")), plot,
                  width = width, height = height, units = "in", device = cairo_pdf)
  if (svg) {
    svg_path <- file.path(FIGURES_DIR, paste0(basename, ".svg"))
    if (requireNamespace("svglite", quietly = TRUE)) {
      # REVISION 2026-10-08: no white background in the SVG. svglite writes a
      # <rect width='100%' height='100%'> that, once the panel is pasted into a larger
      # Inkscape page, becomes an invisible page-sized box covering other panels.
      svg_plot <- if (inherits(plot, "ggplot")) {
        plot + ggplot2::theme(plot.background = ggplot2::element_rect(fill = NA, colour = NA))
      } else plot
      ggplot2::ggsave(svg_path, svg_plot, width = width, height = height, units = "in",
                      bg = "transparent")
      txt <- readLines(svg_path, warn = FALSE)
      txt <- txt[!grepl("<rect width='100%' height='100%'", txt, fixed = TRUE)]
      for (q in c('"', "'")) {
        txt <- gsub(paste0("font-family: ", q, "Liberation Sans", q),
                    paste0("font-family: ", FONT), txt, fixed = TRUE)
      }
      writeLines(txt, svg_path)
    } else {
      warning("svglite not installed; SVG would have text as outlines, skipping")
    }
  }
  invisible(plot)
}

# age- and sex-adjusted Jonckheere-Terpstra
# ordered-trend test on residuals of rank(outcome) ~ Age + Sex.
# Returns list(z, p, n_per_group, n_used).  No external package needed.
jt_test_residuals_age_sex <- function(df, outcome, disease_col, order,
                                       age_col = "Age", sex_col = "Sex") {
  d <- df[complete.cases(df[, c(outcome, disease_col, age_col, sex_col)]), , drop = FALSE]
  d <- d[d[[disease_col]] %in% order, , drop = FALSE]
  if (nrow(d) < 3) return(list(z = NA_real_, p = NA_real_, nper = NULL, n = nrow(d)))
  d$.rk <- rank(d[[outcome]])
  # Residuals of rank ~ Age + Sex; categorical Sex via factor
  fmla <- as.formula(paste0(".rk ~ ", age_col, " + factor(", sex_col, ")"))
  d$.resid <- residuals(lm(fmla, data = d))
  groups <- lapply(order, function(g) d$.resid[d[[disease_col]] == g])
  groups <- groups[lengths(groups) > 0]
  k <- length(groups); N <- sum(lengths(groups))
  U <- 0
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    xi <- groups[[i]]; xj <- groups[[j]]
    U <- U + sum(outer(xj, xi, `>`)) + 0.5 * sum(outer(xj, xi, `==`))
  }
  nis <- lengths(groups)
  EU <- (N^2 - sum(nis^2)) / 4
  VU <- (N^2 * (2*N + 3) - sum(nis^2 * (2*nis + 3))) / 72
  z <- (U - EU) / sqrt(VU)
  p <- 2 * (1 - pnorm(abs(z)))
  list(z = z, p = p,
       nper = setNames(as.integer(nis), order[seq_along(nis)]),
       n = N)
}

# Helper to format the JT label for on-panel annotation.
jt_label <- function(jt) {
  if (is.na(jt$z) || is.na(jt$p)) return(NA_character_)
  p_str <- if (jt$p < 0.001) sprintf("%.1e", jt$p) else sprintf("%.3f", jt$p)
  sprintf("Jonckheere-Terpstra: z=%.2f, p=%s", jt$z, p_str)
}

# Color palette (consistent across all figures)
COLORS <- c(
  "HD" = "#4DBBD5",
  "MGUS" = "#F39B7F",
  "SMM" = "#E64B35",
  "SMM (Untreated)" = "#E64B35",
  "SMM (Treated)" = "#91D1C2",
  "MM" = "#3C5488"
)
