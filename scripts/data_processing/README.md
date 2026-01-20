# Data Processing Scripts

This directory contains scripts for preprocessing the raw data into the formats required for model fitting.

## Files

- **`Prepare_IHME_data.r`**: Processes the IHME COVID-19 data (cases, deaths, vaccination).
- **`Prepare_GISAID_data.r`**: Processes GISAID genomic data to compute variant proportions.
- **`Prepare_travel_data.r`**: Generates the global movement matrix from air travel data.
- **`travel_data_2019.R`**: Prepares travel data specifically for the 2019 baseline (used in future scenario simulations).
- **`Infer_origins.sh`**: Shell script to run the phylogenetic origin inference pipeline.
- **`infer_origins/`**: Subdirectory containing scripts for the phylogenetic inference process.
- **`install_prerequisite.R`**: Helper script to install required R packages.

## Workflow
Typically, these scripts are run in the order described in the main `scripts/README.md` to produce the intermediate data files in `data/`.
