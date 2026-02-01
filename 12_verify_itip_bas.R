# verify_itip_bas.R
# Test script to verify the integration of missForest and BAS in ITIP
# Author: Azman Nads

# setwd("/Users/azmannads/Documents/Research collections/Research extension/ITIP_Research/scripts")

library(testthat)
options(testthat.default_reporter = "summary")
source("04_itip_algorithm.R")

cat("=== ITIP Verification Test ===\n\n")

# 1. Create synthetic data
set.seed(123)
n <- 200
X1 <- rnorm(n)
X2 <- rnorm(n)
# X2 is missing when X1 > 0 (MNAR-ish)
X2_missing <- X2
X2_missing[X1 > 0.5] <- NA

# Outcome depends on X1, X2, and the interaction of X2_imputed * Missing_Indicator
# Let Z = indicator for X2 missing
Z <- as.integer(is.na(X2_missing))
# True model: Y ~ X1 + X2 + 2 * (X2 * Z)
# Note: When X2 is missing, X2_imputed will be close to mean or predicted value
Y_logits <- 0.5 * X1 + 0.5 * X2 + 2.0 * (X2 * Z)
Y_probs <- 1 / (1 + exp(-Y_logits))
Y <- rbinom(n, 1, Y_probs)
X3 <- rnorm(n) # Extra noise variable

data <- data.frame(
    X1 = X1,
    X2 = X2_missing,
    X3 = X3,
    Y = Y
)

cat("Created synthetic dataset (N=200)...\n")

# 2. Run ITIP
cat("Running ITIP algorithm...\n")
# Use a low epsilon to keep interactions, just to test flow
result <- itip(data, outcome = "Y", epsilon = 0.001, verbose = TRUE)

# 3. Verify Outputs
cat("\n=== Verification Results ===\n")

# Check 1: Data Imputed?
if (any(is.na(result$data_imputed))) {
    cat("[FAIL] data_imputed still has missing values.\n")
} else {
    cat("[PASS] missForest imputation successful (no NAs).\n")
}

# Check 2: BAS Results present?
if (is.null(result$bas_results)) {
    cat("[FAIL] BAS results are NULL.\n")
} else {
    cat("[PASS] BAS results generated.\n")
    cat("  Top Features (PIP > 0.5):\n")
    print(result$bas_results$mpm_features)
}

# Check 3: Structure
if ("bas_object" %in% names(result$bas_results)) {
    cat("[PASS] BAS object structure is correct.\n")
} else {
    cat("[FAIL] BAS object missing from results.\n")
}

cat("\n=== End Verification ===\n")
