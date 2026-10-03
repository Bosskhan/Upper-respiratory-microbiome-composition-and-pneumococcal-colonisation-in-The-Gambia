## Script 03: Figure 3 - Dominant genus per participant across visits (alluvial),
##            Figure S3 - treemap of dominant genus across all samples
## Requires: pob_bac_decontam_samples_filtered_asv from script 01, and
##           genus_colours / italic_labels from script 02
##
## Note: this script uses cowplot::get_plot_component()-style legend extraction
## via ggplotGrob(), which is robust across ggplot2 versions including 4.x.
## If using an older ggplot2/cowplot combination, see README for alternatives.

library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(treemapify)
library(cowplot)

pob_final <- pob_bac_decontam_samples_filtered_asv

# ---------------------------------------------------------------------
# 1) Dominant genus per sample
# ---------------------------------------------------------------------

pob_genus_rel <- pob_final %>%
  tax_glom("Genus") %>%
  transform_sample_counts(function(x) x / sum(x))

abund <- as(otu_table(pob_genus_rel), "matrix")
if (taxa_are_rows(pob_genus_rel)) abund <- t(abund)

genus_names <- as.data.frame(tax_table(pob_genus_rel))$Genus
names(genus_names) <- taxa_names(pob_genus_rel)

dominant_asv <- colnames(abund)[max.col(abund, ties.method = "first")]
names(dominant_asv) <- rownames(abund)

dominant_df <- data.frame(sample_data(pob_genus_rel))
dominant_df$dominant_taxa <- unname(genus_names[dominant_asv[rownames(dominant_df)]])
dominant_df <- dominant_df[, c("q3_record_id", "visitmp", "dominant_taxa")]

# Group any genus without a defined colour into "Others"
dominant_df$dominant_taxa <- ifelse(
  dominant_df$dominant_taxa %in% names(genus_colours),
  dominant_df$dominant_taxa,
  "Others"
)

# ---------------------------------------------------------------------
# 2) Alluvial subset: participants with a sample at all 8 visits
# ---------------------------------------------------------------------

alluvial_df <- dominant_df %>%
  group_by(q3_record_id) %>%
  filter(n_distinct(visitmp) == 8) %>%
  ungroup() %>%
  mutate(Year = case_when(
    visitmp %in% paste0("V", 1:4) ~ "Year 1",
    visitmp %in% paste0("V", 5:8) ~ "Year 2"
  ))

genus_levels <- unique(c(
  alluvial_df$dominant_taxa[alluvial_df$Year == "Year 1"],
  alluvial_df$dominant_taxa[alluvial_df$Year == "Year 2"]
))
alluvial_df$dominant_taxa <- factor(alluvial_df$dominant_taxa, levels = genus_levels)

n_participants <- n_distinct(alluvial_df$q3_record_id)

# ---------------------------------------------------------------------
# 3) Panels
# ---------------------------------------------------------------------

plot_alluvial_year <- function(df, title, show_y_axis = TRUE) {
  p <- ggplot(df, aes(x = visitmp, stratum = dominant_taxa,
                      alluvium = q3_record_id, fill = dominant_taxa)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.6) +
    geom_stratum(alpha = 0.8) +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    scale_y_continuous(
      name   = "Proportion of samples",
      breaks = n_participants * c(0, 0.25, 0.5, 0.75, 1),
      labels = function(x) paste0(round(100 * x / n_participants), "%")
    ) +
    scale_fill_manual(values = genus_colours, labels = italic_labels, drop = FALSE) +
    labs(x = "Visit", fill = NULL) +
    ggtitle(title) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title      = element_text(face = "bold", size = 16),
      axis.text       = element_text(face = "bold", size = 14),
      legend.key.size = unit(0.4, "cm"),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 18)
    )

  if (!show_y_axis) {
    p <- p + theme(
      axis.title.y = element_text(colour = "transparent"),
      axis.text.y  = element_text(colour = "transparent"),
      axis.ticks.y = element_blank()
    )
  }
  p
}

fig_year1 <- plot_alluvial_year(filter(alluvial_df, Year == "Year 1"), "Year 1")
fig_year2 <- plot_alluvial_year(filter(alluvial_df, Year == "Year 2"), "Year 2", show_y_axis = FALSE)

# ---------------------------------------------------------------------
# 4) Shared legend (bottom)
# ---------------------------------------------------------------------

legend_df <- data.frame(dominant_taxa = factor(genus_levels, levels = genus_levels))

legend_plot <- ggplot(legend_df, aes(x = 1, y = dominant_taxa, fill = dominant_taxa)) +
  geom_point(shape = 22, size = 6) +
  scale_fill_manual(values = genus_colours, labels = italic_labels, drop = FALSE) +
  guides(fill = guide_legend(nrow = 5, override.aes = list(size = 6))) +
  labs(fill = NULL) +
  theme(legend.text = element_text(size = 13))

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  is_legend <- sapply(g$grobs, function(x) x$name) == "guide-box"
  g$grobs[[which(is_legend)]]
}

shared_legend <- extract_legend(legend_plot + theme(legend.position = "bottom"))

# ---------------------------------------------------------------------
# 5) Combine and save Figure 3
# ---------------------------------------------------------------------

fig3 <- plot_grid(
  plot_grid(fig_year1, fig_year2, ncol = 2, align = "h", axis = "tb"),
  shared_legend,
  ncol = 1, rel_heights = c(1, 0.35)
)

ggsave(
  "output_data/figures/Figure_3_alluvial_dominant_genus.png",
  fig3, width = 16, height = 9, units = "in", dpi = 300, bg = "white"
)

# ---------------------------------------------------------------------
# 6) Figure S3: treemap of dominant genus (all samples)
#    Genera dominant in <1% of samples are grouped as "Others"
# ---------------------------------------------------------------------

treemap_df <- dominant_df %>%
  count(dominant_taxa) %>%
  mutate(group = ifelse(100 * n / sum(n) < 1, "Others", dominant_taxa)) %>%
  group_by(group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    label    = paste0(group, " (", round(100 * n / sum(n)), "%)"),
    fontface = ifelse(group == "Others", "plain", "italic")
  )

fig_treemap <- ggplot(treemap_df, aes(area = n, fill = group, label = label, fontface = fontface)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "centre", grow = FALSE, reflow = TRUE) +
  scale_fill_manual(values = genus_colours) +
  scale_discrete_identity(aesthetics = "fontface", guide = "none") +
  theme(legend.position = "none")

ggsave(
  "output_data/figures/Figure_S3_treemap_dominant_genus.png",
  fig_treemap, width = 14, height = 7, units = "in", dpi = 300, bg = "white"
)
