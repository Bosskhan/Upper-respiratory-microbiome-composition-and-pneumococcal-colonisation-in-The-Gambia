## Script 04: Figure 4 - Longitudinal S. pneumoniae colonisation profiles
## Requires: pob (the original, pre-QC phyloseq object) from script 01,
##           with a 'spn' column (Neg/Pos) in sample_data

library(phyloseq)
library(dplyr)
library(ggplot2)

meta_full <- data.frame(sample_data(pob))
meta_full$visitmp <- factor(meta_full$visitmp, levels = paste0("V", 1:8))

meta_full$Colonisation <- factor(
  meta_full$spn,
  levels = c("Neg", "Pos"),
  labels = c("Non-colonised", "Colonised")
)

carriage_df <- meta_full %>%
  filter(simplifiedtype == "Sample", !is.na(Colonisation)) %>%
  select(q3_record_id, visitmp, Colonisation)

# ---- Participants with a culture result at all 8 visits ----

carriage_8visits <- carriage_df %>%
  group_by(q3_record_id) %>%
  filter(n_distinct(visitmp) == 8) %>%
  ungroup()

# ---- Classify by colonisation frequency ----

carriage_8visits <- carriage_8visits %>%
  group_by(q3_record_id) %>%
  mutate(
    n_colonised = sum(Colonisation == "Colonised"),
    classification = case_when(
      n_colonised %in% 6:8 ~ "High",
      n_colonised %in% 3:5 ~ "Moderate",
      n_colonised %in% 0:2 ~ "Low"
    )
  ) %>%
  ungroup() %>%
  mutate(classification = factor(classification, levels = c("High", "Moderate", "Low")))

# ---- Order participants: by group then by colonisation frequency ----

plot_data <- carriage_8visits %>%
  arrange(classification, desc(n_colonised), q3_record_id) %>%
  mutate(numeric_id = as.numeric(factor(q3_record_id, levels = unique(q3_record_id))))

# ---- Plot ----

fig4 <- ggplot(plot_data, aes(x = visitmp, y = numeric_id, colour = Colonisation)) +
  geom_point(shape = 15, size = 2) +
  facet_grid(classification ~ ., scales = "free_y", space = "free_y") +
  scale_y_continuous(
    trans = "reverse",
    breaks = seq(1, max(plot_data$numeric_id)),
    expand = expansion(add = 0.6)
  ) +
  scale_colour_manual(values = c("Colonised" = "blue", "Non-colonised" = "red")) +
  labs(x = "Visit", y = "Participant", colour = "Colonisation status") +
  theme_bw() +
  theme(
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 11, face = "bold"),
    strip.text.y     = element_text(size = 12, face = "bold", angle = 0),
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

ggsave(
  "output_data/figures/Figure_4_spn_colonisation_profiles.png",
  fig4, width = 10, height = 10, units = "in", dpi = 300, bg = "white"
)
