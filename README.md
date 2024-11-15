#  NCOMMS-24-22914

## Title
Unveiling the pelagic-benthic coupling associated with the biological carbon pump in the Fram Strait (Arctic Ocean)

## Abstract
Settling aggregates transport organic matter from the ocean surface to the deep sea and seafloor. Although plankton communities impact carbon export, the influence of specific organisms and their interactions on export efficiency remains unknown. By examining 15 years of eDNA sequences (18S-V4) from settling and sedimented organic matter in the Fram Strait, we observed that most phylogenetic groups are transferred from pelagic to benthic ecosystems. Notably, organisms such as Chaetoceros socialis, sea-ice diatoms, Radiolaria, and Chaetognatha play key roles in driving the vertical carbon flux to a depth of 200 m, while C. socialis alone is crucial for organic carbon reaching the seafloor. Spatiotemporal changes in community composition indicated that warming events reduce diatom abundance, which can lower the efficiency of this diatom-driven carbon pump. Interestingly, several parasite taxa also displayed strong associations with carbon flux, suggesting that they may induce unknown export mechanisms, potentially affecting pelagic-benthic coupling and ecosystem functioning.
 
## Repository Content

### Purpose
This repository contains R scripts and associated datasets used in the analysis of eDNA sequences for studying the pelagic-benthic coupling and biological carbon pump in the Fram Strait. The code provided here allows for:

Data Import: Reading environmental and omics datasets.
Data Processing: Organizing data matrices for abundance and taxonomy.

1. Data Files
- Pelagic_environmental_data.csv & Benthic_environmental_data.csv: Environmental variables for benthic samples.
- Pelagic_samples.csv & Benthic_samples.csv: Metadata for samples used in the study.
- 18S_Pelagic_data.csv & 18S_Benthic_data.csv: Raw 18S-V4 eDNA sequences and abundance data.

2. Scripts
- Pelagic_ecosystem_analyses & Benthic_ecosystem_analyses.R: Main script for importing, processing, and analyzing data.