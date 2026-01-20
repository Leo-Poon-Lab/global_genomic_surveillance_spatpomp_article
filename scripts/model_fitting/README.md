# Model Fitting Scripts

This directory contains the scripts for calibrating the models to the observed data.

## Files

- **`profiling_HKU/`**: Directory containing the main fitting scripts for the primary SpatPOMP models ($M_0$ and $M_1$). Includes scripts for the 'burn-in' phase (optimization) using high-performance computing resources.
- **`benchmark.R`**: Script for specifying and fitting the benchmark autoregressive statistical model ($M_2$).
- **`build_Omicron20.R`**: Helper script to construct the Omicron model object for fitting.
- **`build_M1.R`**: Helper script to construct the $M_1$ model object for fitting.

## Computation Note
Model fitting for $M_0$ and $M_1$ is computationally intensive and was originally performed on a high-performance computing cluster.
