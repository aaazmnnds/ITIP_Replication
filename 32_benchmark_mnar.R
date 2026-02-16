# 32_benchmark_mnar.R
# Phase 2: Reviewer Response - MNAR Robustness
# Tests ITIP performance when missingness depends on the unobserved value (MNAR)
# Mechanism: logit(P(Z_j=1)) = alpha * X_j + beta

library(MASS)
library(glmnet)
library(missForest)
library(BAS) # Bayesian Adaptive Sampling

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

generate_mnar_data <- function(n = 500, p = 20, rho = 0.5, mnar_strength = 1.5) {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
    X <- mvrnorm(n, rep(0, p), Sigma)
    colnames(X) <- paste0("V", 1:p)

    # Introduce MNAR missingness in V1 and V2
    X_miss <- X

    # MNAR Mechanism: P(Missing V1) depends on V1 itself
    # logit(P) = 1.5 * V1 - 0.5 (as per request)
    miss_prob_v1 <- plogis(mnar_strength * X[, 1] - 0.5)
    missing_idx_v1 <- rbinom(n, 1, miss_prob_v1) == 1
    X_miss[missing_idx_v1, 1] <- NA

    # MNAR Mechanism: P(Missing V2) depends on V2 itself
    miss_prob_v2 <- plogis(mnar_strength * X[, 2] - 0.5)
    missing_idx_v2 <- rbinom(n, 1, miss_prob_v2) == 1
    X_miss[missing_idx_v2, 2] <- NA

    # True Model: Y depends on interaction with missingness
    # 1.5 * (X * (1 - Z)) equivalent to 1.5*X - 1.5*X*Z
    # We want ITIP to find the I term (X*Z)
    lp <- 0.5 * X[, 1] + 0.5 * X[, 2] + 0.3 * X[, 3] +
        1.5 * (X[, 1] * (1 - as.numeric(missing_idx_v1))) +
        1.5 * (X[, 2] * (1 - as.numeric(missing_idx_v2)))

    prob <- plogis(lp)
    y <- rbinom(n, 1, prob)

    return(list(X = X_miss, y = y, true_interactions = c("I_V1", "I_V2")))
}

# --- Main Benchmark Loop ---
p_values <- c(20, 50, 100)
results_list <- list()

for (p in p_values) {
    cat("\nRunning MNAR Benchmark for p =", p, "...\n")

    set.seed(42 + p)
    data <- generate_mnar_data(n = 500, p = p)
    X <- data$X
    y <- data$y
    true_ints <- data$true_interactions

    # 1. ITIP
    cat("  Running ITIP...\n")
    df_itip <- as.data.frame(cbind(X, Y = y))

    # Using n_bins=3 default
    res <- itip(data = df_itip, outcome = "Y", threshold_method = "adaptive", verbose = FALSE)
    itip_ints <- res$pruned_features$kept_interactions

    # 2. LASSO (Baseline) - With Z inclusion fix
    cat("  Running LASSO...\n")
    missing_pattern <- is.na(X)
    imp_lasso <- missForest(X, verbose = FALSE)$ximp
    imp_data <- as.data.frame(imp_lasso)

    # Add Z columns
    miss_cols <- which(colSums(missing_pattern) > 0)
    if (length(miss_cols) > 0) {
        for (k in seq_along(miss_cols)) {
            col_idx <- miss_cols[k]
            col_name <- colnames(X)[col_idx]
            z_name <- paste0("Z_", col_name)
            imp_data[[z_name]] <- as.numeric(missing_pattern[, col_idx])
        }
    }

    # Model Matrix
    X_full <- model.matrix(~ .^2, data = imp_data)[, -1]

    # Add I_extras manually for fair comparison (naming match)
    if (length(miss_cols) > 0) {
        I_extras <- matrix(0, nrow = nrow(X), ncol = length(miss_cols))
        extra_names <- character(length(miss_cols))
        for (k in seq_along(miss_cols)) {
            col_idx <- miss_cols[k]
            col_name <- colnames(X)[col_idx]
            Z_vec <- as.numeric(missing_pattern[, col_idx])
            I_term <- imp_lasso[, col_idx] * Z_vec
            I_extras[, k] <- I_term
            extra_names[k] <- paste0("I_", col_name)
        }
        colnames(I_extras) <- extra_names
        # Combine if not duplicate name issue (simple bind here)
        X_full <- cbind(X_full, I_extras)
    }

    cv_lasso <- cv.glmnet(X_full, y, family = "binomial", alpha = 1)
    coefs <- coef(cv_lasso, s = "lambda.1se")
    selected_lasso <- rownames(coefs)[coefs[, 1] != 0]

    # Clean names
    clean_names <- function(x) {
        if (length(x) == 0) {
            return(character(0))
        }
        x <- gsub("`", "", x)
        is_imp_int <- grepl("^I_", x)
        if (any(!is_imp_int)) {
            x[!is_imp_int] <- sapply(x[!is_imp_int], function(s) {
                parts <- strsplit(s, ":")[[1]]
                paste(sort(parts), collapse = ":")
            })
        }
        return(x)
    }

    lasso_ints <- c(
        grep(":", selected_lasso, value = TRUE),
        grep("^I_", selected_lasso, value = TRUE)
    )

    s_itip <- calc_f1(itip_ints, true_ints)
    s_lasso <- calc_f1(clean_names(lasso_ints), true_ints)

    cat("  Results p=", p, ": ITIP F1=", s_itip["f1"], " LASSO F1=", s_lasso["f1"], "\n")

    results_list[[as.character(p)]] <- c(
        itip_f1 = s_itip["f1"],
        lasso_f1 = s_lasso["f1"]
    )
}



if (dir.exists("results")) {
    saveRDS(results_list, "results/mnar_mnar_results.rds")
} else {
    saveRDS(results_list, "../results/mnar_mnar_results.rds")
}
print(results_list)
