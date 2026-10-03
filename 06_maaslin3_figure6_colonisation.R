## Script 06: MaAsLin3 setup (shared by scripts 06-08) and
##            Figure 6 - MaAsLin3 associations with S. pneumoniae colonisation
## Requires: pob_final (= pob_bac_decontam_samples_filtered_asv) from script 01
##
## Methods: CLR-normalised genus-level abundance/prevalence models,
## participant ID as random effect, read depth included as a covariate,
## age and season controlled for. Significance threshold q < 0.1
## (Benjamini-Hochberg). See manuscript Methods for full description.

library(phyloseq)
library(dplyr)
library(tibble)
library(janitor)
library(maaslin3)
library(ggplot2)

pob_final <- pob_bac_decontam_samples_filtered_asv

# =====================================================================
# MaAsLin3 input setup (shared by scripts 06, 07, 08)
# =====================================================================

pob_genus <- tax_glom(pob_final, taxrank = "Genus")

read_depth_df <- data.frame(
  sampleid = sample_names(pob_genus),
  read_depth = sample_sums(pob_genus)
)

abundance_data <- otu_table(pob_genus)
abundance_data <- t(abundance_data)
abundance_data_df <- as.data.frame(abundance_data)

Genus_names <- as.character(tax_table(pob_genus)[, "Genus"])
colnames(abundance_data_df) <- Genus_names
abundance_data_df <- clean_names(abundance_data_df)

metadata_maaslin <- as(sample_data(pob_genus), "data.frame")
metadata_maaslin <- metadata_maaslin %>%
  rownames_to_column(var = "sample_id") %>%
  left_join(read_depth_df %>% select(sampleid, read_depth), by = c("sample_id" = "sampleid"))

rownames(metadata_maaslin) <- metadata_maaslin$sample_id

# =====================================================================
# Figure 6: S. pneumoniae colonisation
# =====================================================================

metadata_maaslin_spn <- metadata_maaslin %>%
  select(sample_id, spn, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(spn)) %>%
  mutate(
    Colonisation = recode(spn, "Pos" = "Colonised", "Neg" = "Non-colonised"),
    Season       = season_sample_collected
  ) %>%
  select(-spn, -season_sample_collected) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

# Non-colonised set as reference level
metadata_maaslin_spn$Colonisation <- factor(
  metadata_maaslin_spn$Colonisation,
  levels = c("Non-colonised", "Colonised")
)

fit_spn <- maaslin3(
  input_data            = abundance_data_df,
  input_metadata        = metadata_maaslin_spn,
  output                = "MAASLIN_SPN",
  formula               = "~ Colonisation + Season + age_at_visit + read_depth + (1|q3_record_id)",
  normalization         = "CLR",
  transform             = "NONE",
  min_prevalence        = 0.01,
  warn_prevalence       = FALSE,
  augment               = TRUE,
  standardize           = TRUE,
  plot_summary_plot     = TRUE,
  plot_associations     = FALSE,   # avoids a known MaAsLin3/ggplot2 plotting error; model results are unaffected
  max_significance      = 0.1,
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = FALSE,
  max_pngs              = 100,
  summary_plot_first    = 50,
  coef_plot_vars        = c("Colonisation Colonised")
)

# ---- Build Figure 6 from the model results ----

merged_results_spn <- rbind(
  data.frame(fit_spn$fit_data_abundance$results, model_type = "abundance"),
  data.frame(fit_spn$fit_data_prevalence$results, model_type = "prevalence")
)

p6 <- maaslin3:::maaslin3_summary_plot(
  merged_results = merged_results_spn,
  summary_plot_file = "MAASLIN_SPN/figures/summary_plot_custom.png",
  figures_folder     = "MAASLIN_SPN/figures/",
  first_n            = 30,
  max_significance   = 0.1,
  coef_plot_vars     = c("Colonisation Colonised"),
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = FALSE
)

p6$final$scales$scales[[2]]$name <- "Prevalence q-value"
p6$final$scales$scales[[3]]$name <- "Abundance q-value"

fig6 <- p6$final +
  theme(
    axis.text.y      = element_text(face = "italic"),
    strip.text       = element_blank(),
    strip.background = element_blank(),
    plot.margin      = margin(t = 25, r = 5, b = 5, l = 5, unit = "pt")
  ) +
  labs(y = "", x = "MaAsLin3 coefficient") +
  annotate("text", x = -0.6, y = Inf, label = "Non-colonised",
           hjust = 0.5, vjust = -0.5, size = 5, fontface = "bold") +
  annotate("text", x = 0.6, y = Inf, label = "Colonised",
           hjust = 0.5, vjust = -0.5, size = 5, fontface = "bold") +
  coord_cartesian(clip = "off")

ggsave(
  "output_data/figures/Figure_6_maaslin_spn.png",
  fig6, width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)
