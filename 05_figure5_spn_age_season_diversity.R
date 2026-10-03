## Script 05: Figure 5 - S. pneumoniae colonisation by age, season,
##            Shannon diversity, and beta diversity (PCoA)
## Requires: pob and pob_final (= pob_bac_decontam_samples_filtered_asv)
##           from script 01
##
## Note: Figure combination uses patchwork syntax ((p1|p2)/(p3|p4)).
## Developed with ggplot2 4.0.3 / patchwork 1.3.2; see README for notes on
## version compatibility of plot-combination helpers.

library(phyloseq)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)
library(ggsignif)
library(lme4)
library(lmerTest)
library(compositions)
library(vegan)
library(patchwork)

pob_final <- pob_bac_decontam_samples_filtered_asv

# =====================================================================
# Panels a & b: culture-based prevalence by age and season (uses pob,
# the full culture dataset, not pob_final)
# =====================================================================

meta_culture <- data.frame(sample_data(pob)) %>%
  filter(simplifiedtype == "Sample")

# ---- Panel a: age groups ----

meta_age <- meta_culture %>%
  filter(!is.na(age_at_visit)) %>%
  mutate(age_group = cut(
    age_at_visit,
    breaks = c(5, 7, 9, 11, 13, 16),
    labels = c("5-6", "7-8", "9-10", "11-12", "13-15"),
    right = FALSE
  ))

age_summary <- meta_age %>%
  group_by(age_group) %>%
  summarise(n = n(), n_colonised = sum(spn == "Pos"),
            pct = round(100 * n_colonised / n, 1), .groups = "drop")

fig5a <- ggplot(age_summary, aes(x = age_group, y = pct, group = 1)) +
  geom_line(colour = "blue", linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(pct, "%")), vjust = -1, size = 4) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Age in years", y = "Colonisation prevalence") +
  theme_bw(base_size = 16) +
  theme(
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14),
    panel.border = element_blank()
  )

# ---- Panel b: season ----

meta_season <- meta_culture %>%
  filter(!is.na(season_sample_collected), !is.na(spn))

season_summary <- meta_season %>%
  group_by(season_sample_collected) %>%
  summarise(
    n = n(), n_colonised = sum(spn == "Pos"),
    prop = n_colonised / n,
    se   = sqrt(prop * (1 - prop) / n),
    .groups = "drop"
  )

fig5b <- ggplot(season_summary, aes(x = season_sample_collected, y = 100 * prop)) +
  geom_col(fill = "blue", width = 0.2) +
  geom_errorbar(aes(ymin = 100 * (prop - se), ymax = 100 * (prop + se)), width = 0.1) +
  geom_signif(
    comparisons = list(c("Dry season", "Wet season")),
    annotations = "p < 0.001",
    y_position  = 65,
    tip_length  = 0.02,
    textsize    = 6
  ) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Season", y = "Colonisation prevalence") +
  theme_bw(base_size = 16) +
  theme(
    panel.border = element_blank(),
    axis.title   = element_text(face = "bold", size = 16),
    axis.text    = element_text(face = "bold", size = 14)
  )

# ---- Mixed-effects logistic regression: colonisation ~ season + age ----
# Accounts for repeated sampling within participants. p-value for season
# (season_sample_collectedWet season term) is reported in panel b.

meta_model <- meta_culture %>%
  filter(!is.na(age_at_visit), !is.na(season_sample_collected))

meta_model$colonised <- ifelse(meta_model$spn == "Pos", 1, 0)

colonisation_model <- glmer(
  colonised ~ season_sample_collected + age_at_visit + (1 | q3_record_id),
  data = meta_model,
  family = binomial
)

summary(colonisation_model)

# =====================================================================
# Panels c & d: microbiome-based analyses (use pob_final)
# =====================================================================

alpha_div <- estimate_richness(pob_final, measures = "Shannon")

metadata_Shannon <- data.frame(sample_data(pob_final))
metadata_Shannon$Shannon <- alpha_div$Shannon[match(rownames(metadata_Shannon), rownames(alpha_div))]

metadata_Shannon <- metadata_Shannon %>%
  filter(!is.na(spn), !is.na(age_at_visit), !is.na(season_sample_collected))

metadata_Shannon$Colonisation <- factor(
  metadata_Shannon$spn, levels = c("Neg", "Pos"),
  labels = c("Non-colonised", "Colonised")
)

# ---- Panel c: Shannon diversity ----
# Linear mixed-effects model, participant random effect, adjusted for age
# and season. p-value for the spn term is reported in panel c.

spn_model <- lmer(
  Shannon ~ spn + age_at_visit + season_sample_collected + (1 | q3_record_id),
  data = metadata_Shannon
)

summary(spn_model)

fig5c <- ggplot(metadata_Shannon, aes(x = Colonisation, y = Shannon, fill = Colonisation)) +
  geom_boxplot(width = 0.2, outlier.size = 1) +
  geom_signif(
    comparisons = list(c("Non-colonised", "Colonised")),
    annotations = "p < 0.001",
    y_position  = max(metadata_Shannon$Shannon) + 0.2,
    tip_length  = 0.02,
    textsize    = 6
  ) +
  scale_fill_manual(values = c("Non-colonised" = "red", "Colonised" = "blue")) +
  scale_y_continuous(limits = c(2, 7)) +
  labs(x = "S. pneumoniae colonisation", y = "Shannon diversity") +
  theme_bw(base_size = 16) +
  theme(
    panel.border    = element_blank(),
    axis.title      = element_text(face = "bold", size = 16),
    axis.title.x    = element_text(face = "bold.italic", size = 16),
    axis.text       = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# ---- Panel d: Aitchison distance, PCoA, PERMANOVA ----
# Genus-level aggregation; CLR transform with pseudocount of 1, then
# Euclidean distance = Aitchison distance.

pob_genus <- tax_glom(pob_final, taxrank = "Genus")

meta_complete <- data.frame(sample_data(pob_genus)) %>%
  rownames_to_column("sample_id") %>%
  filter(!is.na(age_at_visit), !is.na(season_sample_collected))

pob_sub <- prune_samples(meta_complete$sample_id, pob_genus)

otu_mat <- as(otu_table(pob_sub), "matrix")
if (taxa_are_rows(pob_sub)) otu_mat <- t(otu_mat)

otu_mat_clr <- clr(otu_mat + 1)
dist_aitchison <- dist(otu_mat_clr, method = "euclidean")

meta_complete <- meta_complete[match(rownames(otu_mat), meta_complete$sample_id), ]
stopifnot(identical(meta_complete$sample_id, rownames(otu_mat)))

# PERMANOVA: spn, adjusted for age and season, blocked by participant
set.seed(123)
permanova_res <- adonis2(
  dist_aitchison ~ spn + age_at_visit + season_sample_collected,
  data = meta_complete,
  permutations = 1000,
  by = "margin",
  strata = meta_complete$q3_record_id
)

permanova_res

# PCoA plot
ord_aitch <- ordinate(pob_sub, method = "PCoA", distance = dist_aitchison)

sample_data(pob_sub)$Colonisation <- factor(
  sample_data(pob_sub)$spn, levels = c("Neg", "Pos"),
  labels = c("Non-colonised", "Colonised")
)

fig5d <- plot_ordination(pob_sub, ord_aitch, color = "Colonisation") +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("Non-colonised" = "red", "Colonised" = "blue")) +
  labs(color = bquote(bolditalic("S. pneumoniae") ~ bold("colonisation"))) +
  theme_bw(base_size = 16) +
  theme(
    panel.border    = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 12),
    axis.title      = element_text(face = "bold", size = 16),
    axis.text       = element_text(face = "bold", size = 14)
  )

# =====================================================================
# Combine panels and save
# =====================================================================

p1 <- fig5a + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
p2 <- fig5b + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
p3 <- fig5c + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
p4 <- fig5d + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

combined_spn <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "a")

ggsave(
  "output_data/figures/Figure_5_spn_age_season_diversity.png",
  combined_spn, width = 12, height = 10, dpi = 300, bg = "white"
)
