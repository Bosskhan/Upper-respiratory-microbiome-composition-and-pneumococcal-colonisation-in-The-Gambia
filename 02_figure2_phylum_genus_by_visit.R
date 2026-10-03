## Script 02: Figure 2 - Mean relative abundance by visit, phylum (a) and genus (b)
## Requires: pob_bac_decontam_samples_filtered_asv from script 01

library(phyloseq)
library(dplyr)
library(ggplot2)
library(scales)
library(cowplot)

pob_final <- pob_bac_decontam_samples_filtered_asv

# ---------------------------------------------------------------------
# 1) Colours
# ---------------------------------------------------------------------

phylum_colors <- c(
  Firmicutes       = "#377EB8",
  Proteobacteria   = "#4DAF4A",
  Bacteroidota     = "#880808",
  Fusobacteriota   = "#984EA3",
  Actinobacteriota = "#FF7F00",
  Campilobacterota = "#FFFF33",
  Spirochaetota    = "#A65628",
  Patescibacteria  = "#F781BF"
)

# Genus colours (shades of the phylum colours)
genus_colours <- c(
  "Corynebacterium"             = "#FBCEB1",
  "Actinomyces"                 = "#F2D2BD",
  "Brachybacterium"             = "#FFAC1C",
  "Atopobium"                   = "#CD7F32",
  "f_micrococcaceae"            = "#DAA06D",
  "Cutibacterium"               = "#FF7F00",
  "Moraxella"                   = "#454B1B",
  "Neisseria"                   = "#097969",
  "f_pasteurellaceae"           = "#2E8B57",
  "Haemophilus"                 = "#00FF7F",
  "Pseudomonas"                 = "#E4D00A",
  "Lautropia"                   = "#DFFF00",
  "Actinobacillus"              = "#00FFFF",
  "Aggregatibacter"             = "#50C878",
  "Kingella"                    = "#AFE1AF",
  "f_neisseriaceae"             = "#40E0D0",
  "Campylobacter"               = "#FDDA0D",
  "Streptococcus"               = "#0000FF",
  "Gemella"                     = "#088F8F",
  "Veillonella"                 = "#0096FF",
  "Granulicatella"              = "#7393B3",
  "Megasphaera"                 = "#6F8FAF",
  "Peptostreptococcus"          = "#5D3FD3",
  "Staphylococcus"              = "#00A36C",
  "Lachnoanaerobaculum"         = "#6495ED",
  "Clostridia_UCG-014"          = "#40B5AD",
  "Johnsonella"                 = "#00008B",
  "Oribacterium"                = "#7DF9FF",
  "f_lachnospiraceae"           = "#1434A4",
  "Eubacterium_nodatum_group"   = "#0818A8",
  "Dolosigranulum"              = "#6082B6",
  "Parvimonas"                  = "#0047AB",
  "Catonella"                   = "#3F00FF",
  "f_veillonellaceae"           = "#5F9EA0",
  "Lactococcus"                 = "#ADD8E6",
  "Acholeplasma"                = "#191970",
  "Selenomonas"                 = "#000080",
  "Peptococcus"                 = "#1F51FF",
  "Filifactor"                  = "#A7C7E7",
  "Mogibacterium"               = "#CCCCFF",
  "Stomatobaculum"              = "#B6D0E2",
  "Eubacterium_brachy_group"    = "#96DED1",
  "f_staphylococcaceae"         = "#4169E1",
  "c_clostridia"                = "#0F52BA",
  "Abiotrophia"                 = "#9FE2BF",
  "Aerococcus"                  = "#87CEEB",
  "Fusobacterium"               = "#702963",
  "Leptotrichia"                = "#800080",
  "Streptobacillus"             = "#673147",
  "f_leptotrichiaceae"          = "#AA336A",
  "Absconditabacteriales_SR1"   = "#FAD5A5",
  "Treponema"                   = "#E2DFD2",
  "Prevotella"                  = "#880808",
  "Porphyromonas"               = "#C41E3A",
  "Alloprevotella"              = "#FF0000",
  "Cloacibacterium"             = "#EE4B2B",
  "Capnocytophaga"              = "#A52A2A",
  "Rikenellaceae_RC9_gut_group" = "#800020",
  "Filobacterium"               = "#6E260E",
  "Lentimicrobium"              = "#CC5500",
  "f_prevotellaceae"            = "#E97451",
  "Tannerella"                  = "#D22B2B",
  "Others"                      = "darkgrey"
)

# ---------------------------------------------------------------------
# 2) Helper: mean relative abundance per visit at a given rank
# ---------------------------------------------------------------------

mean_rel_abund_by_visit <- function(ps, rank) {
  ps %>%
    tax_glom(rank) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    merge_samples("visitmp") %>%
    transform_sample_counts(function(x) x / sum(x))
}

theme_fig2 <- theme(
  axis.title      = element_text(size = 14, face = "bold"),
  axis.text.x     = element_text(size = 14, face = "bold", angle = 0, hjust = 0.5),
  axis.text.y     = element_text(size = 14, face = "bold"),
  legend.title    = element_text(size = 16, face = "bold"),
  legend.text     = element_text(size = 14),
  legend.key.size = unit(1, "lines"),
  panel.border    = element_blank()
)

# ---------------------------------------------------------------------
# 3) Panel a: phylum
# ---------------------------------------------------------------------

ps_phylum <- mean_rel_abund_by_visit(pob_final, "Phylum")

fig_phylum <- plot_bar(ps_phylum, fill = "Phylum") +
  scale_fill_manual(values = phylum_colors) +
  scale_y_continuous(labels = percent) +
  labs(x = "Visit", y = "Mean relative abundance", fill = NULL) +
  theme_fig2

# ---------------------------------------------------------------------
# 4) Panel b: genus (top 20, everything else grouped as "Others")
# ---------------------------------------------------------------------

ps_genus <- mean_rel_abund_by_visit(pob_final, "Genus")

genus_df <- psmelt(ps_genus) %>%
  mutate(Genus = as.character(Genus))

top_genera <- genus_df %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  slice_max(total, n = 20) %>%
  pull(Genus)

genus_order <- genus_df %>%
  filter(Genus %in% top_genera) %>%
  group_by(Phylum, Genus) %>%
  summarise(total = sum(Abundance), .groups = "drop") %>%
  arrange(Phylum, desc(total)) %>%
  pull(Genus)

genus_df <- genus_df %>%
  mutate(Genus = factor(
    ifelse(Genus %in% top_genera, Genus, "Others"),
    levels = c(genus_order, "Others")
  ))

# Italic legend labels ("Others" stays plain)
italic_labels <- function(x) {
  parse(text = ifelse(x == "Others", "Others", paste0("italic('", x, "')")))
}

fig_genus <- ggplot(genus_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = genus_colours, labels = italic_labels) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_continuous(labels = percent) +
  labs(x = "Visit", y = "Mean relative abundance", fill = NULL) +
  theme_fig2

# ---------------------------------------------------------------------
# 5) Combine and save
# ---------------------------------------------------------------------

fig2 <- plot_grid(
  fig_phylum, fig_genus,
  labels = c("a", "b"), label_size = 14,
  ncol = 2, align = "h", axis = "tb"
)

ggsave(
  "output_data/figures/Figure_2_genus_phylum_by_visit.png",
  fig2, width = 16, height = 7, units = "in", dpi = 300, bg = "white"
)
