# run_simulation_experiment.R
# Execute Phase 4A: Simulation Study Validation
# Author: Azman Nads

# setwd("/Users/azmannads/Documents/Research collections/Research extension/ITIP_Research/scripts")

library(pROC)
library(dplyr)

source("20_generate_simulation_data.R")
source("04_itip_algorithm.R")

# === Configuration ===
N_REPLICATES <- 10 # Number of simulation runs
SAMPLE_SIZE <- 1000 # Increased sample size for better IG estimation
N_FEATURES <- 10 # Total Features
MISSINGNESS <- 0.3 # 30% Missingness
MECHANISM <- "MNAR" # Difficult case
It_EPSILON <- 0.001 # Lower threshold to detect subtle signal

# Storage for results
results_df <- data.frame(
    run = integer(),
    method = character(),
    auc = numeric(),
    n_features = integer(),
    true_interactions_found = integer(),
    false_interactions_found = integer(),
    stringsAsFactors = FALSE
)

cat(sprintf("=== Starting Simulation Study (%d Replicates) ===\n", N_REPLICATES))
cat(sprintf("Config: N=%d, P=%d, Mechanism=%s, Eps=%.4f\n\n", SAMPLE_SIZE, N_FEATURES, MECHANISM, It_EPSILON))

for (i in 1:N_REPLICATES) {
    set.seed(42 + i)
    if (i %% 1 == 0) cat(sprintf("Replicate %d/%d...\n", i, N_REPLICATES))

    # 1. Generate Data
    sim <- generate_simulation_data(
        n = SAMPLE_SIZE,
        p = N_FEATURES,
        mechanism = MECHANISM,
        seed = 42 + i
    )

    data_train <- sim$data
    true_interactions <- sim$truth$interactions # e.g., "I_X3"

    # 2. Run ITIP
    # We suppress verbose output for the loop (unless debugging needed, keeping verbose=TRUE for now)
    tryCatch(
        {
            itip_res <- itip(data_train, outcome = "Y", epsilon = It_EPSILON, verbose = TRUE)

            # 3. Evaluate Feature selection
            selected_interactions <- itip_res$pruned_features$kept_interactions

            # TPR/FPR
            tp <- sum(selected_interactions %in% true_interactions)
            fp <- sum(!selected_interactions %in% true_interactions)

            # 4. Train Model & Get AUC (In-sample for now, properly should be held-out)
            model <- build_itip_model(itip_res, outcome = "Y", method = "logistic")
            preds <- predict(model, type = "response")
            auc_val <- as.numeric(auc(data_train$Y, preds, quiet = TRUE))

            # Store ITIP results
            results_df <- rbind(results_df, data.frame(
                run = i,
                method = "ITIP",
                auc = auc_val,
                n_features = length(selected_interactions),
                true_interactions_found = tp,
                false_interactions_found = fp
            ))

            # 5. Baseline: Full Model (All interactions)
            tryCatch(
                {
                    if (is.null(itip_res$interactions)) {
                        next
                    }

                    all_interactions <- names(itip_res$interactions)

                    # Base features
                    base_cols <- names(itip_res$data_imputed)[!grepl("^I_", names(itip_res$data_imputed))]
                    base_data <- itip_res$data_imputed[, base_cols, drop = FALSE]

                    # Combine (ensure data frame)
                    interactions_df <- as.data.frame(itip_res$interactions)

                    full_data <- cbind(base_data, interactions_df)
                    full_data$Y <- data_train$Y

                    model_full <- glm(Y ~ ., data = full_data, family = binomial())

                    preds_full <- predict(model_full, type = "response")
                    auc_full <- as.numeric(auc(data_train$Y, preds_full, quiet = TRUE))

                    results_df <- rbind(results_df, data.frame(
                        run = i,
                        method = "Baseline_Full",
                        auc = auc_full,
                        n_features = length(all_interactions),
                        true_interactions_found = 1,
                        false_interactions_found = length(all_interactions) - 1
                    ))
                },
                error = function(e) {
                    cat(sprintf("    Baseline failed in run %d: %s\n", i, e$message))
                }
            )
        },
        error = function(e) {
            cat(sprintf("  Error in run %d: %s\n", i, e$message))
        }
    )
}

# === Summary ===
cat("\n=== Simulation Results Summary ===\n")
summary_stats <- results_df %>%
    group_by(method) %>%
    summarise(
        Mean_AUC = mean(auc),
        SE_AUC = sd(auc) / sqrt(n()),
        Mean_Features = mean(n_features),
        Mean_True_Found = mean(true_interactions_found),
        Mean_False_Found = mean(false_interactions_found)
    )

print(summary_stats)

cat("\nInterpretation:\n")
cat("With lower epsilon and higher N, ITIP should find the true interaction.\n")
