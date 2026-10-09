# ============================================================================
# Purpose:      Figure 5A-B: plasma IL-1B and IL-18 Olink levels, paired pre vs post vaccination, in HD/MGUS/SMM. The single MGUS patient with prior systemic therapy is dropped to keep the cohort treatment-naive. Defines the plot_cytokine() helper reused by Figure5DEF.R.
#               q-values are Benjamini-Hochberg corrected across the whole panel (52 analytes x 3 disease groups = 156 tests), matching the Methods. Effect sizes r are unaffected by the correction.
# Inputs:       data/olink/olink_cytokines.csv.
# Outputs:      figures/Figure5A.png and figures/Figure5B.png (and matching PDFs); side-effect: plot_cytokine() is left in the global environment for Figure5DEF.R.
# Dependencies: R + tidyverse, rstatix, ggpubr; sources ../config.R for FONT and save_figure().
# ============================================================================
source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read Olink proteomics data
data <- read_csv(file.path(DATA_DIR, "olink", "olink_cytokines.csv"))

# drop the single MGUS participant with
# prior systemic therapy from the Olink cohort. The remaining 20/20 SMM and 13/13 MGUS in this
# analysis are treatment-naive (verified by joining the original-ID Olink summary to the ELISA
# master Ever_treated column; mapping kept in scratch, not shared). HD are not applicable.
.OLINK_TREATED_DROP <- c("P49")
data <- data[!data$Patient_ID %in% .OLINK_TREATED_DROP, ]

# Keep only paired samples (pre + post vaccination)
plot_df_index <- data %>%
  group_by(Patient_ID, Cytokine) %>%
  filter(any(Timepoint == "Pre-Vx"), any(Timepoint == "Post-Vx"))

plot_df_index$Timepoint <- ifelse(plot_df_index$Timepoint == "Pre-Vx", "1", "2")
plot_df_index$Timepoint <- factor(plot_df_index$Timepoint, levels = c("1", "2"), ordered = TRUE)
plot_df_index$dummy_timepoint <- as.numeric(plot_df_index$Timepoint)

# Add jitter for paired lines (seeded so the panels are byte-reproducible)
set.seed(2026)
b <- as.list(runif(nrow(plot_df_index), -0.15, 0.15))
b_df <- do.call("rbind", b)
colnames(b_df) <- "add_noise"
plot_df_index <- cbind(plot_df_index, plot_df_index$dummy_timepoint + b_df)
colnames(plot_df_index)[ncol(plot_df_index)] <- "add_noise"

# Ensure paired filtering
master_paired <- plot_df_index %>%
  group_by(Patient_ID) %>%
  filter(any(Timepoint == "1"), any(Timepoint == "2"))

# Build facet labels with sample counts
counts <- as.data.frame(table(master_paired$Timepoint, master_paired$Diagnosis, master_paired$Cytokine))
counts <- counts %>% distinct(Var2, Var3, .keep_all = TRUE)
counts$matchname <- paste0(counts$Var2, "_", counts$Var3)
master_paired$matchname <- paste0(master_paired$Diagnosis, "_", master_paired$Cytokine)
# Display label: "Healthy" is shown as HD throughout the figure.
# REVISION 2026-10-08: facet titles are plotmath expressions so that n is italic.
counts$n_ind <- paste0("atop('", ifelse(counts$Var2 == "Healthy", "HD", as.character(counts$Var2)),
                       "', '('*italic(n)*'=", counts$Freq, ")')")
master_paired$n_ind <- counts$n_ind[match(master_paired$matchname, counts$matchname)]

# Panel-wide multiple-testing family.
# The Methods state that paired pre- versus post-vaccination comparisons are BH-corrected
# "across all cytokines tested", so the correction is computed ONCE over the whole panel
# rather than within each analyte. The family is every analyte with at least 3 paired
# participants in all three disease groups, which resolves to the same 52 analytes carried
# into the Figure 2D cross-sectional screen, giving 52 x 3 = 156 tests. Restricting to
# analytes with enough paired data is required because a signed-rank test on fewer than 3
# pairs is uninformative and would only inflate the family.
.panel_analytes <- master_paired %>%
  group_by(Cytokine, Diagnosis) %>%
  summarise(n_pt = n_distinct(Patient_ID), .groups = "drop") %>%
  group_by(Cytokine) %>%
  filter(n() == 3, all(n_pt >= 3)) %>%
  pull(Cytokine) %>%
  unique()

PANEL_Q <- master_paired %>%
  filter(Cytokine %in% .panel_analytes) %>%
  group_by(Cytokine, Diagnosis) %>%
  wilcox_test(Level ~ Timepoint, paired = TRUE) %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  select(Cytokine, Diagnosis, p.adj)
stopifnot(length(.panel_analytes) == 52, nrow(PANEL_Q) == 156)
message(sprintf("panel-wide BH family: %d analytes x 3 groups = %d tests",
                length(.panel_analytes), nrow(PANEL_Q)))

# Display names for y-axis labels
cyto_labels <- c("IL1B" = "IL-1\u03B2", "IL18" = "IL-18",
                 "NUB1" = "NUB1", "DDX58" = "DDX58", "MMP7" = "MMP7")

# Helper function for individual cytokine plots ---
plot_cytokine <- function(cyto_name, paired_data) {
  plot_paired_df <- paired_data %>% filter(Cytokine %in% cyto_name)
  plot_paired_df$VaccineTimepoint <- plot_paired_df$Timepoint
  y_label <- ifelse(cyto_name %in% names(cyto_labels), cyto_labels[cyto_name], cyto_name)

  p <- ggplot(plot_paired_df, aes(x = VaccineTimepoint, y = Level, fill = VaccineTimepoint)) +
    geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, show.legend = FALSE, color = NA) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(add_noise, y = Level, group = Patient_ID),
              color = "black", size = 0.3, alpha = 0.5, show.legend = FALSE) +
    geom_pointrange(aes(add_noise, ymin = Level, ymax = Level),
                    shape = 21, size = 0.5, alpha = 1, show.legend = FALSE, linetype = "dotted") +
    scale_fill_manual(values = c("steelblue", "tomato2")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, color = "black"),
          # REVISION 2026-10-08: larger text (was 16 / 18 / 18), identical in Figure 5C.
          axis.text = element_text(size = 18, color = "black", family = FONT),
          # one-line y-axis titles: 15 pt is the largest size at which the longest title
          # ("DDX58 plasma levels (Olink)") fits the axis height
          axis.title = element_text(size = 15, color = "black", family = FONT),
          strip.text = element_text(size = 20, face = "plain", family = FONT),
          axis.title.x = element_blank(),
          strip.background = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_x_discrete(breaks = c(1, 2), labels = c("Pre", "Post")) +
    facet_wrap(. ~ n_ind, labeller = label_parsed) +
    # REVISION 2026-10-08: values are Olink plasma levels (linear scale), not NPX.
    ylab(paste0(y_label, " plasma levels (Olink)")) + xlab("")

  # Paired Wilcoxon signed-rank test per disease group. The q-value is NOT corrected
  # here: it is looked up from PANEL_Q, the panel-wide BH family of 156 tests.
  stat.test <- plot_paired_df %>%
    group_by(n_ind) %>%
    wilcox_test(Level ~ VaccineTimepoint, paired = TRUE) %>%
    add_significance()

  stat.test <- stat.test %>%
    left_join(plot_paired_df %>% ungroup() %>% distinct(n_ind, Diagnosis), by = "n_ind") %>%
    left_join(PANEL_Q %>% filter(Cytokine == cyto_name) %>% select(Diagnosis, p.adj),
              by = "Diagnosis")
  stopifnot(!any(is.na(stat.test$p.adj)))

  # Paired rank-biserial effect sizes per disease group
  eff <- plot_paired_df %>%
    group_by(n_ind) %>%
    wilcox_effsize(Level ~ VaccineTimepoint, paired = TRUE)

  stat.test <- stat.test %>%
    add_xy_position(x = "VaccineTimepoint") %>%
    left_join(eff %>% select(n_ind, effsize), by = "n_ind") %>%
    mutate(
      q_num = signif(p.adj, 2),
      # REVISION 2026-10-08: plotmath labels so that q and r are italic, numbers roman.
      q_fmt = paste0("italic(q)*'=", q_num, "'"),
      bracket_label = ifelse(q_num < 0.1,
                             paste0(q_fmt, "*', '*italic(r)*'=", round(abs(effsize), 2), "'"),
                             q_fmt))

  # JT subtitle removed per user instruction; the figure now shows
  # only the per-disease paired Wilcoxon q (with effect size r when q<0.1). Bracket label
  # text color forced to black (was gray under ggpubr default).
  final <- p +
    # REVISION 2026-10-08: geom_bracket() directly, because stat_pvalue_manual() does not pass
    # type = "expression" through and the plotmath labels would print as raw text.
    geom_bracket(data = stat.test, aes(xmin = xmin, xmax = xmax, y.position = y.position,
                                       label = bracket_label),
                 type = "expression", family = FONT, tip.length = 0.02, bracket.nudge.y = 1,
                 label.size = 4.6, size = 0.3, color = "black", inherit.aes = FALSE) +
    coord_cartesian(clip = "off") +
    theme(panel.spacing.x = unit(0.4, "lines"),
          plot.margin = margin(8, 12, 4, 4)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)))

  return(final)
}

# Figure 5A: IL-1B
p_il1b <- plot_cytokine("IL1B", master_paired)
save_figure(p_il1b, "Figure5A", width = 6, height = 4)

# Figure 5B: IL-18
p_il18 <- plot_cytokine("IL18", master_paired)
save_figure(p_il18, "Figure5B", width = 6, height = 4)
