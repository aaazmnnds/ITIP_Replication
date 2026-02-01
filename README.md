# ITIP Replication Package

This repository contains the R code to reproduce the results and figures for the paper:
**"Information-Theoretic Interaction Pruning for Machine Learning with Missing Data"**

## Directory Structure

*   `00_setup_environment.R`: Installs necessary dependencies.
*   `03_bas_utils.R`, `04_itip_algorithm.R`: Core algorithm implementation.
*   `14_validate_prop4.R`: Reproduces the empirical validation of the BAS proxy (Appendix B, Figure 3).
*   `22_reproduce_toy_example.R`: Reproduces the "Illustrative Example" (Section 5, Table 2, Figure 2).
*   `30_benchmark_baselines.R`: Reproduces the large-scale simulation benchmarks (Section 7, Table 3).
*   `41_real_world_validation.R`: Code for the MIMIC-III application (Section 8).

## Prerequisites

All scripts are written in R (v4.3+). Install dependencies by running:
```R
source("00_setup_environment.R")
```
Key packages: `missForest`, `glmnet`, `stabs`, `BAS`, `MASS`.

## Reproduction Instructions

### 1. Illustrative Example (Section 5)
To generate the conditional information gain table (Table 2) and the pruning behavior plot (Figure 2):
```bash
Rscript 22_reproduce_toy_example.R
```
*   **Expected Output:** A table of IG scores for 5 interactions; Plot saved to `../results/`.
*   **Note:** Uses fixed seed 3 for exact reproducibility.

### 2. Simulation Benchmarks (Section 7)
To generate the performance comparison table (Table 3) for p=20, 50, 100:
```bash
Rscript 30_benchmark_baselines.R
```
*   **Expected Output:** F1, Precision, and Recall scores for ITIP, LASSO, and Stability Selection across dimensions.
*   **Runtime:** Approximately 1-2 hours depending on hardware (due to iterative simulations).
*   **Parameters:** n=500, rho=0.5, beta=0.8, epsilon=0.005.

### 3. Proposition Validation (Appendix)
To validate the correlation between BAS-proxy and Mutual Information:
```bash
Rscript 14_validate_prop4.R
```

## Data Availability
The simulation data is generated synthetically within the scripts.
The MIMIC-III data used in Section 8 requires a signed data use agreement and cannot be shared. The script `41_real_world_validation.R` provides the analysis logic assuming the data is locally available in `../data/`. Users must obtain access via [PhysioNet](https://physionet.org/content/mimiciii/) to run this specific script.

## Author & Contact

**Azman Nads**  
Informatics and Data Science Program  
Graduate School of Advanced Science and Engineering  
Hiroshima University, Japan  
Email: `azmannads@msutawi-tawi.edu.ph`

For questions regarding the code or the paper, please open a GitHub issue or contact the author via email.
