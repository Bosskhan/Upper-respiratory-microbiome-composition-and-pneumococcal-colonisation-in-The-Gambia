## Script 08: MaAsLin3 models for additional covariates (Table S3)
## Requires: abundance_data_df and metadata_maaslin from script 06
##
## Each covariate is a participant-/household-level survey item recorded
## once per participant rather than at every visit; tidyr::fill() is used
## to propagate the recorded value to all of that participant's visit-rows
## before fitting. All models share the same formula structure (covariate +
## season + age + read depth + participant random effect) and the same
## min_prevalence/warn_prevalence settings, matching the Methods description.

library(dplyr)
library(tidyr)
library(tibble)
library(maaslin3)

run_maaslin_covariate <- function(metadata, covariate_formula_term, output_dir,
                                   max_pngs = 100, summary_plot_first = 50) {
  maaslin3(
    input_data            = abundance_data_df,
    input_metadata        = metadata,
    output                = output_dir,
    formula               = paste0("~ ", covariate_formula_term,
                                   " + season_sample_collected + age_at_visit + read_depth + (1|q3_record_id)"),
    normalization         = "CLR",
    transform             = "NONE",
    min_prevalence        = 0.01,
    warn_prevalence       = FALSE,
    augment               = TRUE,
    standardize           = TRUE,
    plot_summary_plot     = TRUE,
    max_significance      = 0.1,
    median_comparison_abundance  = TRUE,
    median_comparison_prevalence = FALSE,
    max_pngs              = max_pngs,
    summary_plot_first    = summary_plot_first
  )
}

# ---- Waste disposal ----

metadata_maaslin_waste <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(q_waste_disposal, .direction = "downup") %>%
  ungroup() %>%
  mutate(q_waste_disposal = recode(q_waste_disposal,
    "Open burning in a designated area in thecommunity" = "open burning outside the household")) %>%
  select(sample_id, q_waste_disposal, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(q_waste_disposal != "Other", !is.na(q_waste_disposal)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_waste <- run_maaslin_covariate(metadata_maaslin_waste, "q_waste_disposal", "MAASLIN_WASTE")

# ---- PCV vaccination status ----

metadata_maaslin_pcv <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(q19_pcv, .direction = "downup") %>%
  ungroup() %>%
  select(sample_id, q19_pcv, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(q19_pcv)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_pcv <- run_maaslin_covariate(metadata_maaslin_pcv, "q19_pcv", "MAASLIN_PCV")

# ---- Fuel source (Charcoal vs Firewood only; households reporting
#      both or neither fuel type are excluded from this analysis) ----

metadata_maaslin_fuel <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(q1_flsrc_1, q1_flsrc_2, .direction = "downup") %>%
  ungroup() %>%
  mutate(fuel_source = case_when(
    q1_flsrc_1 == 1 & q1_flsrc_2 == 0 ~ "Charcoal",
    q1_flsrc_1 == 0 & q1_flsrc_2 == 1 ~ "Firewood",
    TRUE ~ NA_character_
  )) %>%
  select(sample_id, fuel_source, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(fuel_source)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_fuel <- run_maaslin_covariate(metadata_maaslin_fuel, "fuel_source", "MAASLIN_FUEL")

# ---- Pets (dog and/or cat in the household) ----

metadata_maaslin_pets <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(a_dg, q2b_ct, .direction = "downup") %>%
  ungroup() %>%
  mutate(Pets = case_when(
    a_dg == 1 & q2b_ct == 2 ~ "Yes",
    a_dg == 2 & q2b_ct == 1 ~ "Yes",
    a_dg == 1 & q2b_ct == 1 ~ "Yes",
    a_dg == 2 & q2b_ct == 2 ~ "No"
  )) %>%
  select(sample_id, Pets, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(Pets)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_pets <- run_maaslin_covariate(metadata_maaslin_pets, "Pets", "MAASLIN_PETS")

# ---- Livestock (cattle, goats, or sheep in the household) ----

metadata_maaslin_livestock <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(q2d_cw, q2e_gt, q2f_shp, .direction = "downup") %>%
  ungroup() %>%
  mutate(Livestock = if_else(q2d_cw == 1 | q2e_gt == 1 | q2f_shp == 1, "Yes", "No")) %>%
  select(sample_id, Livestock, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(Livestock)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_livestock <- run_maaslin_covariate(metadata_maaslin_livestock, "Livestock", "MAASLIN_LIVESTOCK")

# ---- Sex ----

metadata_maaslin_sex <- metadata_maaslin %>%
  group_by(q3_record_id) %>%
  fill(q12_sex, .direction = "downup") %>%
  ungroup() %>%
  select(sample_id, q12_sex, season_sample_collected, age_at_visit, read_depth, q3_record_id) %>%
  filter(!is.na(q12_sex)) %>%
  as.data.frame() %>%
  column_to_rownames("sample_id")

fit_sex <- run_maaslin_covariate(metadata_maaslin_sex, "q12_sex", "MAASLIN_SEX")

# Significant results for each covariate (q < 0.1) are written automatically
# by maaslin3 to <output_dir>/significant_results.tsv and compiled into
# Table S3 for the manuscript.
