# 34_sensitivity_bins.R
# Phase 2: Reviewer Response - Discretization Sensitivity
# Tests ITIP stability across different bin counts (n_bins = 3, 5, 10)

library(missForest)
library(BAS)

# Source ITIP functions
if (file.exists("03_bas_utils.R")) {
    source("03_bas_utils.R")
    source("04_itip_algorithm.R")
} else {
    source("scripts/03_bas_utils.R")
    source("scripts/04_itip_algorithm.R")
}

# --- Helper ---
calc_jaccard <- function(set1, set2) {
    if (length(set1) == 0 && length(set2) == 0) {
        return(1)
    }
    return(length(intersect(set1, set2)) / length(union(set1, set2)))
}

# --- Data Generation (Reusing standard simulation) ---
generate_data_sensitivity <- function(n = 500, p = 20) {
    # Simple independent features for speed
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("V", 1:p)

    # Missingness
    miss_prob <- plogis(X[, 3] - 0.5)
    missing_idx <- rbinom(n, 1, miss_prob) == 1
    X[missing_idx, 1] <- NA

    # True Model: Y ~ V1 + V2 + 1.5*(V1*Z1)
    # n_bins shouldn't affect V1*Z1 detection heavily ideally
    Z1 <- as.numeric(missing_idx)
    lp <- 0.5 * X[, 1] + 0.5 * X[, 2] + 1.5 * (X[, 1] * (1 - Z1)) # 1-Z1 implies outcome depends on observed val
    # Actually ITIP looks for X_imp * Z.
    # Let's use the explicit structure: Y depends on V1 when missing.
    # lp = ... + 1.5 * V1 * Z1

    # Correct mechanism for ITIP target:
    # Interaction I_V1 = V1_imp * Z1
    X_imp <- X
    X_imp[is.na(X)] <- 0 # dummy
    lp <- 0.5 * X[, 2] + 1.5 * (X_imp[, 1] * Z1)

    y <- rbinom(n, 1, plogis(lp))

    return(list(X = X, y = y))
}

# --- Sensitivity Loop ---
bin_settings <- c(3, 5, 10)
results <- list()

cat("Running Discretization Sensitivity Check...\n")
set.seed(123)
data <- generate_data_sensitivity(n = 500, p = 20)
df <- as.data.frame(data$X)
df$Y <- data$y

sets <- list()

for (b in bin_settings) {
    cat("  Testing n_bins =", b, "...\n")
    # Capture output to keep clean
    capture.output({
        res <- itip(df, "Y", threshold_method = "adaptive", n_bins = b, verbose = FALSE)
    })

    selected <- res$pruned_features$kept_interactions
    sets[[as.character(b)]] <- selected
    cat("    Selected:", length(selected), "interactions\n")
    print(selected)
}

# Compare consistency
cat("\nConsistency (Jaccard Similarity):\n")
j_3_5 <- calc_jaccard(sets[["3"]], sets[["5"]])
j_5_10 <- calc_jaccard(sets[["5"]], sets[["10"]])
j_3_10 <- calc_jaccard(sets[["3"]], sets[["10"]])

cat("  bins 3 vs 5: ", j_3_5, "\n")
cat("  bins 5 vs 10:", j_5_10, "\n")
cat("  bins 3 vs 10:", j_3_10, "\n")

results <- list(sets = sets, jaccard = c(j35 = j_3_5, j510 = j_5_10, j310 = j_3_10))

if (dir.exists("results")) {
    saveRDS(results, "results/sensitivity_bins_results.rds")
} else {
    saveRDS(results, "../results/sensitivity_bins_results.rds")
}
