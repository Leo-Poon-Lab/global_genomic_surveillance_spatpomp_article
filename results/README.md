# Results Documentation

This directory contains the outputs from the model calibration, simulation, and analysis.

## Directory Structure

### `model_data/`
Contains intermediate data files and processed model outputs.
- **Content**: 
    - Calibrated model objects.
    - Simulated trajectories for different scenarios.
    - Bootstrapped results for strategy evaluation.
- **Usage**: These files are loaded by the analysis scripts in `scripts/final_analysis` to reproduce the figures and tables in the manuscript.

### `benchmark/`
Contains results from the benchmark models.
- **Content**: Comparison metrics (likelihoods, AIC, BIC) between the primary mechanistic model ($M_0$), simplified mechanistic model ($M_1$), and statistical model ($M_2$).
