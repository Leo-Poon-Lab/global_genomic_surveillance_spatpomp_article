# Model Building Scripts

This directory contains the R scripts used to define the structure of the SpatPOMP models used in the study.

## Files

- **`spatPomp_Omicron20.R`**: Defines the primary multi-strain metapopulation model ($M_0$). This model includes phylogenetic origin data in the observation process.
- **`spatPomp_M1.R`**: Defines the simplified mechanistic model ($M_1$) which excludes phylogenetically-inferred sequence origin data. Used for validation.
- **`spatPomp_M3.R`**: Defines the retrospective intervention model ($M_3$), an extension of $M_0$ for simulating alternative surveillance strategies during the Omicron wave.
- **`spatPomp_M4.R`**: Defines the future scenario model ($M_4$), an extension of $M_0$ for simulating hypothetical future variants under 2019 travel conditions.
- **`global_data_prepare.R`**: Script for preparing global data structures required by the models.

## Usage
These functions are sourced by the fitting and simulation scripts in `../model_fitting` and `../model_simulation`.
