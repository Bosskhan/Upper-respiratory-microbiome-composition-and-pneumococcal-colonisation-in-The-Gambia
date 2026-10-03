## S. pneumoniae and the upper respiratory tract microbiome: 16S rRNA analysis
## Script 01: Build phyloseq object, taxa filtering, contaminant removal,
##            and genus name cleaning
##
## Input:  data/asv_table.xlsx, data/taxonomy_table.xlsx, data/metadata.xlsx
##         (metadata.xlsx is not included in this repository - see README)
## Output: pob_bac_decontam_samples_filtered_asv (used by all later scripts)

library(writexl)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(decontam)
library(phyloseq)

# =====================================================================
# SECTION A: Building phyloseq object, taxa filtering, contaminant removal
# =====================================================================

# ---- 1) Building phyloseq ----

asv_raw <- read_excel("data/asv_table.xlsx")
asv <- as.data.frame(asv_raw)
rownames(asv) <- asv$ASV
asv$ASV <- NULL
asv <- as.matrix(asv)

tax_raw <- read_excel("data/taxonomy_table.xlsx")
tax <- as.data.frame(tax_raw)
rownames(tax) <- tax$ASV
tax$ASV <- NULL
tax <- as.matrix(tax)

metadata_raw <- read_excel("data/metadata.xlsx")
metadata <- as.data.frame(metadata_raw)
rownames(metadata) <- metadata$sequencing_id
metadata$sequencing_id <- NULL

pob <- phyloseq(
  otu_table(asv, taxa_are_rows = TRUE),
  tax_table(tax),
  sample_data(metadata)
)

# ---- 2) Subset bacteria domain, remove other taxa at Phylum level and
#         Chloroplast contaminant from the data ----

pob_bac <- subset_taxa(
  pob,
  Taxon == "Bacteria" &
    !is.na(Phylum) & Phylum != "" & Phylum != "NA" &
    (is.na(Order) | Order != "Chloroplast")
)

# Rename ASVs for easier reporting
asv_ids <- taxa_names(pob_bac)

asv_id_map <- data.frame(
  Original_ASV = asv_ids,
  ASV = paste0("ASV", seq_along(asv_ids)),
  stringsAsFactors = FALSE
)
taxa_names(pob_bac) <- asv_id_map$ASV

write_xlsx(asv_id_map, "output_data/spreadsheet/asv_id_map.xlsx")

# ---- 3) Identify contaminant ASVs using decontam ----

sample_data(pob_bac)$is.neg <-
  sample_data(pob_bac)$simplifiedtype == "NC"

contamdf.prev <- isContaminant(
  pob_bac,
  method = "prevalence",
  neg = "is.neg",
  threshold = 0.10
)

contaminant_asv <- rownames(contamdf.prev)[contamdf.prev$contaminant]

contaminant_taxa <- as.data.frame(tax_table(pob_bac)[contaminant_asv, ])
contaminant_taxa$ASV <- rownames(contaminant_taxa)
contaminant_taxa <- contaminant_taxa[c("ASV", setdiff(names(contaminant_taxa), "ASV"))]

write_xlsx(contaminant_taxa, "output_data/spreadsheet/decontam_contaminant_ASVs.xlsx")

# Remove contaminant ASVs
pob_bac_decontam <- prune_taxa(
  !taxa_names(pob_bac) %in% contaminant_asv,
  pob_bac
)

# ---- Figure S1: Library sizes after decontam contaminant removal ----

df <- data.frame(sample_data(pob_bac_decontam))
df$LibrarySize <- sample_sums(pob_bac_decontam)

df_samples <- df[as.character(df$simplifiedtype) == "Sample", , drop = FALSE]
df_controls <- df[
  as.character(df$type) %in% c("Extraction negative control", "Field control", "PCR NTC"),
  , drop = FALSE
]

df_samples <- df_samples[order(df_samples$LibrarySize), ]
df_controls <- df_controls[order(df_controls$LibrarySize), ]
df_samples$Index <- seq_len(nrow(df_samples))
df_controls$Index <- seq_len(nrow(df_controls))

df_controls$type <- factor(
  df_controls$type,
  levels = c("Extraction negative control", "Field control", "PCR NTC"),
  labels = c("DNA Extraction negative control", "Field control", "PCR no template control")
)

plot_samples <- ggplot(df_samples, aes(x = Index, y = LibrarySize)) +
  geom_point(color = "steelblue") +
  labs(title = "Samples", x = "Samples", y = "Read count") +
  scale_y_continuous(breaks = seq(0, 14000, by = 2000)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

plot_controls <- ggplot(df_controls, aes(x = Index, y = LibrarySize, color = type)) +
  geom_point(size = 3) +
  labs(title = "Controls", x = "Controls", y = NULL) +
  scale_y_continuous(breaks = seq(0, 1000, by = 200)) +
  scale_color_manual(values = c(
    "DNA Extraction negative control" = "red",
    "Field control" = "orange",
    "PCR no template control" = "navy"
  )) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank()
  )

figure_S1 <- plot_grid(plot_samples, plot_controls, ncol = 2, rel_widths = c(2, 1), align = "h")

ggsave(
  "output_data/figures/Figure_S1_library_size_controls.png",
  figure_S1, width = 14, height = 7, units = "in", dpi = 300
)

# Retain biological samples only
pob_bac_decontam_samples <- prune_samples(
  sample_data(pob_bac_decontam)$simplifiedtype == "Sample",
  pob_bac_decontam
)

# ---- 4) Remove samples with fewer than 1000 reads ----

sample_reads <- sample_sums(pob_bac_decontam_samples)

pob_bac_decontam_samples_filtered <- prune_samples(
  sample_reads >= 1000,
  pob_bac_decontam_samples
)

# ---- 5) Filter ASVs by abundance and prevalence ----
# Keep ASVs reaching >=0.1% relative abundance in >=1% of samples

pob_bac_decontam_samples_filtered_rel <- transform_sample_counts(
  pob_bac_decontam_samples_filtered,
  function(x) x / sum(x)
)

asv_prevalence <- apply(
  otu_table(pob_bac_decontam_samples_filtered_rel), 1,
  function(x) sum(x >= 0.001)
)

asvs_to_keep <- names(asv_prevalence)[
  asv_prevalence >= ceiling(0.01 * nsamples(pob_bac_decontam_samples_filtered))
]

pob_bac_decontam_samples_filtered_asv <- prune_taxa(
  asvs_to_keep,
  pob_bac_decontam_samples_filtered
)

# ---- Figure S2: samples analysed in the final dataset per participant ----

meta_all <- data.frame(sample_data(pob))
meta_all$sample_id <- rownames(meta_all)
meta_all <- meta_all[meta_all$simplifiedtype %in% "Sample", ]

meta_all$Status <- ifelse(
  meta_all$sample_id %in% sample_names(pob_bac_decontam_samples_filtered_asv),
  "Analysed", "Not analysed"
)

meta_all$visitmp <- factor(meta_all$visitmp)

n_analysed <- tapply(meta_all$Status == "Analysed", meta_all$q3_record_id, sum)
participant_order <- names(n_analysed)[order(-n_analysed, names(n_analysed))]

participant_key <- data.frame(
  q3_record_id = participant_order,
  participant_number = paste0("P", seq_along(participant_order)),
  stringsAsFactors = FALSE
)

write_xlsx(participant_key, "output_data/spreadsheet/participant_id_key.xlsx")

meta_all$participant_number <- participant_key$participant_number[
  match(meta_all$q3_record_id, participant_key$q3_record_id)
]

plot_df <- expand.grid(
  participant_number = participant_key$participant_number,
  visitmp = levels(meta_all$visitmp),
  stringsAsFactors = FALSE
)

plot_df$Status <- meta_all$Status[match(
  paste(plot_df$participant_number, plot_df$visitmp),
  paste(meta_all$participant_number, meta_all$visitmp)
)]
plot_df$Status[is.na(plot_df$Status)] <- "Not collected"

plot_df$Status <- factor(plot_df$Status, levels = c("Analysed", "Not analysed", "Not collected"))
plot_df$visitmp <- factor(plot_df$visitmp, levels = levels(meta_all$visitmp))
plot_df$participant_number <- factor(plot_df$participant_number, levels = rev(participant_key$participant_number))

status_colours <- c(
  "Not collected" = "#00A896",
  "Not analysed"  = "#FF5A5F",
  "Analysed"      = "#FFB400"
)

figure_S2 <- ggplot(plot_df, aes(x = visitmp, y = participant_number, fill = Status)) +
  geom_point(shape = 21, colour = "grey55", size = 2.5, stroke = 0.3) +
  scale_fill_manual(values = status_colours, drop = FALSE) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  labs(x = "Visit", y = "Participant", fill = NULL) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 11)
  )

ggsave(
  "output_data/figures/Figure_S2_sample_inclusion_by_visit.png",
  figure_S2, width = 9, height = 17, units = "in", dpi = 300
)

# ---- 6) Clean genus names ----

# a. Remove ASV classified as Mitochondria
tax_check <- as.data.frame(tax_table(pob_bac_decontam_samples_filtered_asv))
mito_asv <- rownames(tax_check)[tax_check$Family %in% "Mitochondria"]

pob_bac_decontam_samples_filtered_asv <- prune_taxa(
  !taxa_names(pob_bac_decontam_samples_filtered_asv) %in% mito_asv,
  pob_bac_decontam_samples_filtered_asv
)

# b. Clean placeholder genus names (e.g. "Prevotellaceae Genus" -> "f_prevotellaceae")
tax_df <- as.data.frame(tax_table(pob_bac_decontam_samples_filtered_asv), stringsAsFactors = FALSE)
genus_original <- tax_df$Genus

needs_fix <- grepl(" Genus$", tax_df$Genus) | grepl("aceae$", tax_df$Genus)
parent <- sub(" Genus$", "", tax_df$Genus)

same <- function(x, y) !is.na(y) & x == y
prefix <- ifelse(same(parent, tax_df$Family), "f_",
          ifelse(same(parent, tax_df$Order),  "o_",
          ifelse(same(parent, tax_df$Class),  "c_",
          ifelse(same(parent, tax_df$Phylum), "p_", "?_"))))

tax_df$Genus[needs_fix] <- tolower(paste0(prefix, parent))[needs_fix]
stopifnot(!any(grepl("^\\?_", tax_df$Genus)))

genus_map <- data.frame(
  ASV = rownames(tax_df),
  Genus_original = genus_original,
  Genus_clean = tax_df$Genus,
  stringsAsFactors = FALSE
)

# Strip brackets/parentheses from polyphyletic group names
tax_df$Genus <- gsub("[][()]", "", tax_df$Genus)

tax_table(pob_bac_decontam_samples_filtered_asv) <- tax_table(as.matrix(tax_df))

write_xlsx(genus_map, "output_data/spreadsheet/genus_name_map.xlsx")

# `pob` and `pob_bac_decontam_samples_filtered_asv` are used by all subsequent scripts.
# Save them if running scripts in separate sessions:
# saveRDS(pob, "output_data/pob.rds")
# saveRDS(pob_bac_decontam_samples_filtered_asv, "output_data/pob_final.rds")
