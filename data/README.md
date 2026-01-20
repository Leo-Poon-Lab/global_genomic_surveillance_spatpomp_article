# Data Documentation

This directory contains the data sources used for the project's analysis and modeling.

## Directory Structure

### `IHME`
Contains daily COVID-19 data sourced from the Institute for Health Metrics and Evaluation (IHME).
- **Contents**: Daily reported cases, reported deaths (with excess mortality scalars), infection detection ratios (IDR), infection fatality ratios (IFR), and vaccination coverage.
- **Usage**: Used for parameterizing the epidemiological model and calibrating surveillance efforts.

### `gisaid_data`
Contains SARS-CoV-2 genomic data.
- **Source**: GISAID (Global Initiative on Sharing All Influenza Data).
- **Contents**: Daily sequenced cases for various lineages (Delta, Omicron BA.1, Omicron BA.2, Others).
- **Usage**: Used to calculate the Diagnostic-to-Sequencing Ratio (DSR) and track variant proportions.
- **Note**: GISAID data is subject to their user agreement and may not be fully redistributable. Readers should refer to the original source for access (Identifier: EPI_SET_240923dt).

### `open_sky_data`, `openflights`, `our_airports_data`
Contains air travel and airport data.
- **Usage**: Combined to generate the high-resolution global passenger travel matrix.

### `covid-policy-tracker`
Contains international travel control measures.
- **Source**: Oxford COVID-19 Government Response Tracker.
- **Usage**: Accounts for policy-driven changes in travel and transmission dynamics.

### `estimated_movement_matrix`
- **Description**: The generated global movement matrix between regions based on the air travel datasets.

### `transit_data`
- **Description**: Data relating to transit passenger volumes at major airports.

### `inferred_origins`
- **Description**: Results of the phylogenetic analysis, inferring the origin of sequenced cases.
- **Usage**: Integrated into the observation model to constrain transmission dynamics.

## Data Availability
Some raw data (e.g., specific GISAID sequences or proprietary flight data) may not be publicly available in this repository due to licensing restrictions. Please refer to the manuscript's Data Availability section and the original sources for access.
