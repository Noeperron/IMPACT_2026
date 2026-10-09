# ============================================================================
# Purpose:      Supplementary Figure 5: sample shipping (shipped vs not-shipped)
#               does not confound SARS-CoV-2 spike-specific TCR clonotype
#               frequencies in treatment-naive SMM, in either pre- or
#               post-vaccination samples. Paired-design control panel for
#               Figure 4B.
# Inputs:       data/tcr/tcr_clonotype_proportions.rds (de-identified ClusTCR
#               per-patient clonotype proportions, bootstrap iterations)
#               data/tcr/tcr_shipping_status.csv (de-identified shipping flag
#               per patient x timepoint, SMM-untreated only)
# Outputs:      figures/SupFig5.png (and matching PDF + SVG)
# Dependencies: R + dplyr, tidyr, readr, ggplot2, ggpubr, rstatix, plotrix;
#               sources ../config.R
# ============================================================================

source("../config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(plotrix)
})

# Load de-identified inputs ---------------------------------------------
clonotype_prop_df <- readRDS(file.path(DATA_DIR, "tcr", "tcr_clonotype_proportions.rds"))
shipping_df <- read_csv(file.path(DATA_DIR, "tcr", "tcr_shipping_status.csv"),
                        show_col_types = FALSE)

# Build the paired-per-patient mean + SD COVID-clonotype dataframe ------
# Mirrors the plot_paired_df construction in the legacy Figure4A-C script:
# aggregate bootstrap iterations to per-patient mean and SD per timepoint.
average_prop <- clonotype_prop_df %>%
  group_by(PatientID, VaccineTimepoint, DiseaseStatus) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)), .groups = "drop")

# Keep CEF / COVID summaries
average_prop_covid <- average_prop %>%
  select(PatientID, VaccineTimepoint, DiseaseStatus,
         pre_COVID_CT_prop_mean, post_COVID_CT_prop_mean,
         pre_COVID_CT_prop_sd, post_COVID_CT_prop_sd)

# COVID mean + SD per (PatientID, VaccineTimepoint). VaccineTimepoint 1 uses the
# pre_* columns; VaccineTimepoint 2 uses the post_* columns.
plot_paired_df <- average_prop_covid %>%
  mutate(
    COVID_mean_prop = ifelse(VaccineTimepoint == 1, pre_COVID_CT_prop_mean, post_COVID_CT_prop_mean),
    COVID_SD_prop   = ifelse(VaccineTimepoint == 1, pre_COVID_CT_prop_sd,   post_COVID_CT_prop_sd)
  ) %>%
  select(PatientID, DiseaseStatus, VaccineTimepoint, COVID_mean_prop, COVID_SD_prop)

# Express clonotype proportions as percentages so this panel is on the same scale
# as Figure 4B, whose y-axis is also a percentage (01_Figure4B.R applies the same
# x100 conversion). Statistics are unaffected by the constant scale factor.
plot_paired_df <- plot_paired_df %>%
  mutate(across(c(COVID_mean_prop, COVID_SD_prop), ~ . * 100))

# Join with shipping-status table ---------------------------------------
shipping_df <- shipping_df %>%
  mutate(VaccineTimepoint = as.integer(VaccineTimepoint),
         Shipping_Status = as.factor(Shipping_Status))

smm_ship_paired <- plot_paired_df %>%
  filter(DiseaseStatus == "SMM:Untreated") %>%
  inner_join(shipping_df, by = c("PatientID", "DiseaseStatus", "VaccineTimepoint")) %>%
  mutate(VaccineTimepoint = as.factor(VaccineTimepoint))

# Statistics: Wilcoxon (shipped vs not) per timepoint, BH adjusted ------
stat_ship <- smm_ship_paired %>%
  group_by(VaccineTimepoint) %>%
  wilcox_test(COVID_mean_prop ~ Shipping_Status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_y_position(formula = COVID_mean_prop ~ Shipping_Status,
                 data = smm_ship_paired, fun = "max", step.increase = 0.1)

# n labels, counted separately for each timepoint. These previously came from the
# pre-vaccination facet alone and were then reused for both facets via scale_x_discrete,
# which mislabelled the post-vaccination not-shipped group as n=5 when it is n=6.
shipment_counts <- smm_ship_paired %>%
  count(VaccineTimepoint, Shipping_Status) %>%
  # REVISION 2026-10-08: plotmath label so that n is italic.
  mutate(n_label = paste0("italic(n)*'=", n, "'"))

timepoint_labels <- c("1" = "Pre-Vx", "2" = "Post-Vx")

# Plot ------------------------------------------------------------------
p_ship <- smm_ship_paired %>%
  ggplot(aes(x = Shipping_Status, y = COVID_mean_prop, fill = Shipping_Status)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.75),
              linewidth = 1, show.legend = FALSE, color = NA) +
  geom_boxplot(lwd = 1, alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_pointrange(aes(ymin = COVID_mean_prop - COVID_SD_prop,
                      ymax = COVID_mean_prop + COVID_SD_prop),
                  linetype = "dotted",
                  position = position_jitter(width = 0.2),
                  shape = 21, show.legend = FALSE) +
  scale_fill_manual(values = c("steelblue", "tomato2")) +
  scale_x_discrete(labels = c("1" = "Shipped", "2" = "Not shipped")) +
  # drawn as a layer rather than baked into the axis labels, so each facet shows its own n
  geom_text(data = shipment_counts,
            aes(x = Shipping_Status, y = -Inf, label = n_label),
            inherit.aes = FALSE, vjust = -0.9, size = 3.9, family = FONT, parse = TRUE) +
  facet_grid(~ VaccineTimepoint, scales = "free",
             labeller = labeller(VaccineTimepoint = timepoint_labels)) +
  # REVISION 2026-10-08: geom_bracket() with a plotmath label so that p is italic;
  # stat_pvalue_manual() does not pass type = "expression" through (see 01_Figure5AB.R).
  geom_bracket(data = stat_ship %>% mutate(p_label = paste0("italic(p)*' = ", p, "'")),
               aes(xmin = group1, xmax = group2, y.position = y.position, label = p_label),
               type = "expression", family = FONT, tip.length = 0.03, label.size = 3.88,
               inherit.aes = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 11, family = FONT),
        axis.text.y = element_text(color = "black", size = 12, family = FONT),
        axis.title = element_text(size = 14, family = FONT),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "plain", size = 14, family = FONT),
        plot.title = element_text(size = 14, hjust = 0.5, family = FONT),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 10, unit = "pt")) +
  xlab("") +
  ylab("SARS-CoV-2 spike-specific clonotypes (%)") +
  ggtitle("SMM (treatment-naive): shipped vs not-shipped")

# Save ------------------------------------------------------------------
save_figure(p_ship, "SupFig5", width = 7, height = 4.5)
print(stat_ship)
