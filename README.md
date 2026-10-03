# Gambia-URT-microbiome-pneumococcal-colonisation

Upper respiratory tract bacterial microbiome features associated with
*Streptococcus pneumoniae* colonisation and environmental factors among
healthy older children in a rural setting in The Gambia.

**Status:** Manuscript submitted to BMC Microbiology.

## Summary

This project characterises the upper respiratory tract (URT) bacterial
microbiome of healthy children aged 5–14 years in The Gambia, sampled
longitudinally over eight visits, using 16S rRNA sequencing alongside
culture-based detection of *S. pneumoniae*. The analysis explores
microbiome composition and temporal dynamics in this age group, and
their relationship to pneumococcal colonisation status, age, and season.

## Repository contents

This repository contains analysis code only. 
The sequence data generated in this study have been deposited in the
SRA under the BioProject accession PRJNA1538550.The dataset used for
the microbiome composition and dominant genus analyses in this paper
is available on Figshare at https://doi.org/10.6084/m9.figshare.34056576.
The metadata are available upon request, given the
requirement for participant confidentiality under the informed consent
obtained for this study. Access can be requested through the Gambia
Government/MRC Joint Ethics Committee, facilitated by MRC Unit The
Gambia (http://www.mrc.gm/) through Martin Antonio
(manotonio@mrc.gm) or Dam Khan (dam.khan@lshtm.ac.uk).

- `scripts/` — R scripts, numbered in the order they should be run,
  covering phyloseq construction and quality control, all main and
  supplementary figures, and other analyses
- `session_info.txt` — R and package versions used to generate the
  results in the manuscript.

### Script order

| Script | Produces |
|---|---|
| `01_phyloseq_setup_and_filtering.R` | phyloseq object construction, decontam contaminant removal, read-depth and prevalence filtering, genus name cleaning, Figures S1–S2 |
| `02_figure2_phylum_genus_by_visit.R` | Figure 2 |
| `03_figure3_alluvial_and_figureS3_treemap.R` | Figure 3, Figure S3 |
| `04_figure4_spn_colonisation_profiles.R` | Figure 4 |
| `05_figure5_spn_age_season_diversity.R` | Figure 5 |
| `06_maaslin3_figure6_colonisation.R` | MaAsLin3 setup (shared by scripts 06–08), Figure 6 |
| `07_maaslin3_figure7_season_age.R` | Figure 7 |
| `08_maaslin3_covariates_table_S3.R` | Table S3 (additional covariate associations) |

Script `01` must be run first; it builds the filtered phyloseq object
(`pob_bac_decontam_samples_filtered_asv`) and the unfiltered object
(`pob`) used by all later scripts. Script `06` builds shared MaAsLin3
input objects used by scripts `07` and `08`. Scripts within each of
these two groups should be run in the same R session, or their output
objects saved and reloaded between sessions.

## Requirements

R version 4.4.2 was used for the manuscript analysis. Required
packages are loaded at the top of each script; see `session_info.txt`
for exact versions used. Key packages: `phyloseq`, `decontam`,
`maaslin3`, `lme4`/`lmerTest`, `vegan`, `compositions`, `ggplot2`,
`cowplot`, `patchwork`, `ggalluvial`, `treemapify`.




Questions about this project can be directed to dam.khan@lshtm.ac.uk.
