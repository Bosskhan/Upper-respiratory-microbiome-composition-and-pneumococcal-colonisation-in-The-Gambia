## Script 07: Figure 7 - MaAsLin3 associations with age and season
## Requires: abundance_data_df and metadata_maaslin from script 06

library(dplyr)
library(tibble)
library(maaslin3)
library(ggplot2)

metadata_maaslin_season_age <- metadata_maaslin %>%
  select(sample_id, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(season_sample_collected), !is.na(age_at_visit)) %>%
  mutate(
    Season = case_when(
      season_sample_collected == "Dry season" ~ "Dry",
      season_sample_collected == "Wet season" ~ "Wet",
      TRUE ~ season_sample_collected
    ),
    Age = age_at_visit
  ) %>%
  select(-season_sample_collected, -age_at_visit) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_season_age <- maaslin3(
  input_data            = abundance_data_df,
  input_metadata        = metadata_maaslin_season_age,
  output                = "MAASLIN_SEASON_AGE",
  formula               = "~ Age + Season + read_depth + (1|q3_record_id)",
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
  summary_plot_first    = 30,
  heatmap_vars          = c("Age", "Season Wet")
)

# ---- Build Figure 7 from the model results ----

merged_results_season_age <- rbind(
  data.frame(fit_season_age$fit_data_abundance$results, model_type = "abundance"),
  data.frame(fit_season_age$fit_data_prevalence$results, model_type = "prevalence")
)

p7 <- maaslin3:::maaslin3_summary_plot(
  merged_results      = merged_results_season_age,
  summary_plot_file   = "MAASLIN_SEASON_AGE/figures/summary_plot_custom.png",
  figures_folder      = "MAASLIN_SEASON_AGE/figures/",
  first_n             = 30,
  max_significance    = 0.1,
  heatmap_vars        = c("Age", "Season Wet"),
  median_comparison_abundance  = TRUE,
  median_comparison_prevalence = FALSE
)

fig7 <- p7$final +
  labs(y = NULL) +
  guides(
    fill   = guide_legend(title = "MaAsLin3 coefficient"),
    colour = guide_legend(title = "q-value")
  ) +
  scale_x_discrete(labels = c("Age" = "Age", "Season Wet" = "Wet Season")) +
  theme(axis.text.y = element_text(face = "italic", size = 12))

ggsave(
  "output_data/figures/Figure_7_maaslin_season_age.png",
  fig7, width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)
