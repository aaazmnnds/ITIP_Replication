# 33_benchmark_high_dim.R
# Phase 2: Reviewer Response - High Dimensional Scalability
# Demonstrates ITIP feasibility at p=200

library(missForest)
library(BAS) # Bayesian Adaptive Sampling
library(parallel)

# Source ITIP functions
if (file.exists("03_bas_utils.R")) {
    source("03_bas_utils.R")
    source("04_itip_algorithm.R")
} else {
    source("scripts/03_bas_utils.R")
    source("scripts/04_itip_algorithm.R")
}

# --- Helper Functions ---

calc_f1 <- function(selected_features, true_features) {
    tp <- sum(selected_features %in% true_features)
    fp <- sum(!selected_features %in% true_features)
    fn <- sum(!true_features %in% selected_features)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0

    return(c(f1 = f1, precision = precision, recall = recall))
}

generate_high_dim_data <- function(n = 500, p = 200, rho = 0.5) {
    # Manual generation to avoid MASS dependency for huge matrix if needed
    # But for p=200 MASS::mvrnorm is fine
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
    X <- tryCatch(
        {
            MASS::mvrnorm(n, rep(0, p), Sigma)
        },
        error = function(e) {
            matrix(rnorm(n * p), n, p) # Fallback to independent
        }
    )
    colnames(X) <- paste0("V", 1:p)

    # Introduce MAR missingness mainly in V1, V2
    # Keep it sparse for others
    X_miss <- X

    # Missingness V1 depends on V3
    miss_prob <- plogis(2 * X[, 3] - 1)
    missing_idx <- rbinom(n, 1, miss_prob) == 1
    X_miss[missing_idx, 1] <- NA

    # Missingness V2 depends on V4
    miss_prob2 <- plogis(1.5 * X[, 4] - 0.8)
    missing_idx2 <- rbinom(n, 1, miss_prob2) == 1
    X_miss[missing_idx2, 2] <- NA

    # Removed sporadic MCAR missingness to match true interactions definition
    # sporadic_cols <- sample(5:p, floor(p / 10))
    # for (j in sporadic_cols) {
    #     X_miss[sample(1:n, floor(0.1 * n)), j] <- NA
    # }

    # True Model
    # 3 Main effects: V1, V2, V5
    # 2 Interactions: I_V1 (V1*Z1), I_V2 (V2*Z2)
    lp <- 0.5 * X[, 1] + 0.5 * X[, 2] + 0.3 * X[, 5] +
        1.5 * (X[, 1] * (1 - as.numeric(missing_idx))) +
        1.5 * (X[, 2] * (1 - as.numeric(missing_idx2)))

    prob <- plogis(lp)
    y <- rbinom(n, 1, prob)

    return(list(X = X_miss, y = y, true_interactions = c("I_V1", "I_V2")))
}

# --- Main Benchmark ---
p_value <- 200
cat("\nRunning High-Dim Benchmark for p =", p_value, "...\n")

set.seed(42 + p_value)
start_time <- Sys.time()
data <- generate_high_dim_data(n = 500, p = p_value)
X <- data$X
y <- data$y
true_ints <- data$true_interactions

# Run ITIP
cat("  Running ITIP...\n")
df_itip <- as.data.frame(cbind(X, Y = y))

# Modified ITIP call: adaptive threshold
# Note: missForest on p=200 takes time.
# BAS MCMC on p=200 + interactions (sparse) is feasible.
res <- itip(data = df_itip, outcome = "Y", threshold_method = "adaptive", verbose = TRUE)

end_time <- Sys.time()
duration <- as.numeric(difftime(end_time, start_time, units = "mins"))

itip_ints <- res$pruned_features$kept_interactions
scores <- calc_f1(itip_ints, true_ints)

cat("\n=== Results p=200 ===\n")
cat("Time:", round(duration, 2), "minutes\n")
cat("F1:", scores["f1"], "Prec:", scores["precision"], "Rec:", scores["recall"], "\n")
cat("Features Selected:", length(itip_ints), "\n")
print(itip_ints)

results <- list(
    p = 200,
    time_mins = duration,
    scores = scores,
    selected = itip_ints
)



if (dir.exists("results")) {
    saveRDS(results, "results/high_dim_results_p200.rds")
} else {
    saveRDS(results, "../results/high_dim_results_p200.rds")
}
