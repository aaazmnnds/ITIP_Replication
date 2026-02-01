library(parallel)
library(BAS)
library(missForest)
library(glmnet)

source("20_generate_simulation_data.R")

# Extended Data Generation to support Sparse Large-Scale Simulations
generate_large_scale_data <- function(n, p, n_true_int = 3, miss_rate = 0.3, mechanism = "MAR") {
    set.seed(runif(1) * 10000)

    # Features
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- paste0("X", seq_len(p))

    # Missingness
    Z <- matrix(0, n, p)
    for (j in seq_len(p)) {
        if (mechanism == "MCAR") {
            Z[, j] <- rbinom(n, 1, miss_rate)
        } else {
            # MAR: Logic for MAR
            prob_miss <- plogis(-1.5 + 0.5 * X[, (j %% p) + 1]) # depends on neighbor
            Z[, j] <- rbinom(n, 1, prob_miss)
        }
    }

    # Outcome Y
    # 3 Main effects, 3 Interactions
    main_idx <- 1:3
    int_idx <- 4:6

    # Interaction: X_j * Z_j
    logits <- 0.5 * X[, main_idx[1]] - 0.5 * X[, main_idx[2]] + 0.3 * X[, main_idx[3]]
    for (idx in int_idx) {
        logits <- logits + 1.5 * (X[, idx] * Z[, idx])
    }

    Y <- rbinom(n, 1, plogis(logits))

    X_obs <- X
    X_obs[Z == 1] <- NA
    X_df <- as.data.frame(X_obs)
    X_df$Y <- Y

    return(list(data = X_df, Z = Z, true_interactions = paste0("X", int_idx)))
}

# Metrics
calculate_f1 <- function(selected, true_set) {
    tp <- sum(selected %in% true_set)
    fp <- sum(!(selected %in% true_set))
    fn <- sum(!(true_set %in% selected))

    precision <- if (tp + fp > 0) tp / (tp + fp) else 0
    recall <- if (tp + fn > 0) tp / (tp + fn) else 0

    f1 <- if (precision + recall > 0) 2 * (precision * recall) / (precision + recall) else 0
    return(f1)
}

# Run single iteration
run_iteration <- function(n, p) {
    sim <- generate_large_scale_data(n, p)

    # 1. Impute (Faster missForest)
    imp <- missForest(sim$data, ntree = 10)$ximp

    # 2. Construct Interactions (X_imp * Z)
    int_data <- imp
    for (j in seq_len(p)) {
        col_name <- paste0("X", j)
        if (any(sim$Z[, j] == 1)) {
            int_col <- imp[[col_name]] * sim$Z[, j]
            int_data[[paste0(col_name, "_Z")]] <- int_col
        }
    }

    # 3. ITIP (BAS) - Faster BAS
    formula_full <- as.formula(paste("Y ~ ."))
    model_bas <- bas.glm(formula_full, data = int_data, family = binomial(), method = "BAS", update = 200)

    # Threshold selection (Prop 4 proxy)
    selected_itip <- names(int_data)[model_bas$probne0 > 0.5]
    # We only care about interactions for this benchmark
    selected_int_itip <- selected_itip[grepl("_Z", selected_itip)]
    selected_int_itip <- gsub("_Z", "", selected_int_itip)

    # 4. Compare with LASSO (on imputed data WITH interactions)
    X_mat <- model.matrix(Y ~ . - 1, data = int_data)
    cv_lasso <- cv.glmnet(X_mat, sim$data$Y, family = "binomial")
    lasso_coef <- as.matrix(coef(cv_lasso, s = "lambda.min"))
    selected_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
    selected_int_lasso <- selected_lasso[grepl("_Z", selected_lasso)]
    selected_int_lasso <- gsub("_Z", "", selected_int_lasso)

    # Calculate F1
    f1_itip <- calculate_f1(selected_int_itip, sim$true_interactions)
    f1_lasso <- calculate_f1(selected_int_lasso, sim$true_interactions)

    return(c(itip = f1_itip, lasso = f1_lasso))
}

# Run experiment
p_values <- c(20, 50, 100)
results_summary <- list()

for (p in p_values) {
    cat("Running p =", p, "\n")
    res <- replicate(3, run_iteration(n = 500, p = p))
    results_summary[[as.character(p)]] <- rowMeans(res)
}

print(results_summary)
