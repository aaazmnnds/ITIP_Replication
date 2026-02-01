# ITIP Research - Proposition 4 Validation Script

## Objective
The JMLR reviewer flagged Proposition 4 (using $-\log(1-\gamma_j)$ as an IG proxy) as a heuristic. This script will:
1. Simulate a set of variables with varying degrees of association (mutual information) with the outcome.
2. Estimate the Information Gain (IG) using a standard non-parametric estimator (e.g., KNN-based CMI).
3. Compute the BAS posterior inclusion probability $(\gamma_j)$ and its proxy value $-\log(1-\gamma_j)$.
4. Analyze the correlation and monotonicity between the proxy and the true IG.

## Setup
- Use `BAS` package for posterior probabilities.
- Use `infotheo` or `entropy` package for "true" Mutual Information.
- Generate data with varied signal-to-noise ratios (SNR).

## Expected Output
- A plot showing $\IG$ on the x-axis and Proxy on the y-axis.
- Spearman correlation coefficient.
