# comprehensive_benchmarks.R
# Final JMLR-ready Horse-Race: ITIP vs. LASSO vs. Elastic Net vs. Stability Selection vs. Hierarchical LASSO
# Author: ITIP Research Team

library(BAS)
library(missForest)
library(glmnet)
library(MASS) # For mvrnorm

# Source new baselines
source("30_benchmark_baselines.R")

# Helper to normalize feature names for comparison
normalize_name <- function(x) {
    x <- gsub("[^[:alnum:]]", "", x) # Keep only alphanumeric
    return(tolower(x))
}

# Benchmarking Function
run_benchmark <- function(n = 500, p = 20, n_true_int = 3, miss_rate = 0.3, debug = FALSE) {
    set.seed(NULL)

    # 1. Generate Correlated Features (rho=0.3 to match manuscript)
    rho <- 0.3
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    colnames(X) <- paste0("X", seq_len(p))

    # 2. Missingness
    Z <- matrix(rbinom(n * p, 1, miss_rate), n, p)
    colnames(Z) <- paste0("Z", seq_len(p))

    # Signal: Interactions X1:Z1, X2:Z2, X3:Z3 + Main Effects X1, X2
    int_idx <- 1:n_true_int
    logits <- 1.0 * X[, 1] - 1.0 * X[, 2]
    for (idx in int_idx) {
        logits <- logits + 5.0 * (X[, idx] * Z[, idx]) # Strong interaction for detection
    }
    Y <- as.numeric(rbinom(n, 1, plogis(logits)))

    X_obs <- X
    X_obs[Z == 1] <- NA
    X_df <- as.data.frame(X_obs)

    # 3. Imputation
    imp_res <- missForest(X_df, ntree = 20, verbose = FALSE)
    imp <- imp_res$ximp

    # 4. Feature Matrix Construction
    X_int <- imp
    for (j in seq_len(p)) {
        X_int[[paste0("Z", j)]] <- Z[, j]
        X_int[[paste0("X", j, "_Z", j)]] <- imp[[paste0("X", j)]] * Z[, j]
    }

    # Matrix for regularized baselines (Saturated Model)
    X_mat <- model.matrix(~ . - 1, data = X_int)
    X_base <- cbind(imp, Z)

    # --- METHODS ---

    # 1. ITIP (Algorithm 1: Pruning)
    selected_int_itip <- c()
    itip_pips <- numeric(p)
    for (j in seq_len(p)) {
        # Minimal conditional model: Y ~ Zj + Ij
        df_j <- data.frame(Y = Y, Zj = Z[, j], Ij = X_int[[paste0("X", j, "_Z", j)]])
        # include.index=1 forces Zj (intercept is always excluded from probne0)
        fit_j <- tryCatch(
            {
                bas.glm(Y ~ Zj + Ij,
                    data = df_j, family = binomial(),
                    include.index = c(1), n.models = 4, method = "BAS", update = 100
                )
            },
            error = function(e) {
                return(NULL)
            }
        )

        if (!is.null(fit_j)) {
            itip_pips[j] <- fit_j$probne0[3] # PIP for Ij
            if (itip_pips[j] > 0.5) {
                selected_int_itip <- c(selected_int_itip, paste0("X", j, "_Z", j))
            }
        }
    }

    # 2. LASSO (Saturated)
    cv_lasso <- cv.glmnet(X_mat, Y, family = "binomial", alpha = 1)
    lasso_coef <- as.matrix(coef(cv_lasso, s = "lambda.min"))[-1, 1]
    selected_lasso <- names(lasso_coef)[lasso_coef != 0]
    selected_int_lasso <- selected_lasso[grepl("_Z", selected_lasso)]

    # 3. Elastic Net (Saturated)
    cv_en <- cv.glmnet(X_mat, Y, family = "binomial", alpha = 0.5)
    en_coef <- as.matrix(coef(cv_en, s = "lambda.min"))[-1, 1]
    selected_en <- names(en_coef)[en_coef != 0]
    selected_int_en <- selected_en[grepl("_Z", selected_en)]

    # 4. Stability Selection (Saturated)
    selected_stabs <- run_stability_selection(X_mat, Y)
    selected_int_stabs <- selected_stabs[grepl("_Z", selected_stabs)]

    # 5. Hierarchical LASSO
    selected_int_hier <- run_hierarchical_lasso(X_base, Y)
    # Filter to only keep X:Z or Z:X interactions
    selected_int_hier <- selected_int_hier[grepl("X.*:Z|Z.*:X", selected_int_hier)]

    # --- METRICS ---
    true_set_norm <- normalize_name(paste0("X", int_idx, "Z", int_idx))

    calc_f1_norm <- function(sel_names, tru_norm) {
        # Find all terms that look like interactions
        sel_ints <- sel_names[grepl("[_:]", sel_names)]
        sel_norm <- normalize_name(sel_ints)

        tp <- sum(sel_norm %in% tru_norm)
        fp <- sum(!(sel_norm %in% tru_norm))
        fn <- sum(!(tru_norm %in% sel_norm))

        prec <- if (tp + fp > 0) tp / (tp + fp) else 0
        rec <- if (tp + fn > 0) tp / (tp + fn) else 0
        return(if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0)
    }

    if (debug) {
        cat("\nDEBUG INFO (p =", p, "):\n")
        cat("Y Balance:", round(mean(Y), 3), "\n")
        cat("True Interactions:", paste(true_set_norm, collapse = ", "), "\n")
        cat("ITIP Selected:", paste(selected_int_itip, collapse = ", "), "\n")
        cat("LASSO Selected:", paste(selected_int_lasso, collapse = ", "), "\n")
        cat("hierNet Selected:", paste(selected_int_hier, collapse = ", "), "\n")
    }

    return(c(
        itip = calc_f1_norm(selected_int_itip, true_set_norm),
        lasso = calc_f1_norm(selected_int_lasso, true_set_norm),
        en = calc_f1_norm(selected_int_en, true_set_norm),
        stabs = calc_f1_norm(selected_int_stabs, true_set_norm),
        hierNet = calc_f1_norm(selected_int_hier, true_set_norm)
    ))
}

# Run Experiment
if (!interactive()) {
    # Ensure results directory
    if (!dir.exists("../results")) dir.create("../results")

    p_values <- c(20, 50, 100)
    n_reps <- 10
    final_results <- list()

    for (p in p_values) {
        cat(sprintf("\n--- Benchmarking p = %d (n = 500, reps = %d) ---\n", p, n_reps))
        res_list <- list()
        for (i in 1:n_reps) {
            cat(sprintf("  Rep %d/%d...\n", i, n_reps))
            res_list[[i]] <- run_benchmark(p = p, debug = (i == 1))
        }
        res_df <- do.call(rbind, res_list)
        final_results[[as.character(p)]] <- colMeans(res_df)
        saveRDS(final_results, "../results/comprehensive_benchmarks.rds")
        print(round(final_results[[as.character(p)]], 3))
    }

    cat("\nFinal Results Summary:\n")
    print(final_results)
}
