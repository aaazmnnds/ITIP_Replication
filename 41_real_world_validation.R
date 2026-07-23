# 41_real_world_validation.R
# Leakage-safe, per-fold Stage 1/Stage 2 ITIP evaluation with paired testing on MIMIC-III
# This script evaluates the complete two-stage ITIP framework (Entropy Screening + Nested LRT)
# inside a strict 5-fold cross-validation loop to avoid information leakage.

library(glmnet)
library(pROC)
library(missForest)
source("04_itip_algorithm.R")

cat("Loading MIMIC Septic Shock Data...\n")
data <- read.csv("../data/mimic_septic_shock_final.csv")
data <- data[, !(names(data) %in% c("HADM_ID"))]

missing_rates <- colMeans(is.na(data))
keep_cols <- names(missing_rates[missing_rates < 1])
data <- data[, keep_cols]

y <- data$mortality
X_miss <- data[, !(names(data) %in% c("mortality"))]

# Create fixed CV folds
set.seed(42)
n <- nrow(X_miss)
folds <- sample(rep(1:5, length.out = n))

# Results storage
metrics_std <- list(auc = numeric(5), f1 = numeric(5), brier = numeric(5))
metrics_full <- list(auc = numeric(5), f1 = numeric(5), brier = numeric(5))
metrics_s1 <- list(auc = numeric(5), f1 = numeric(5), brier = numeric(5))
metrics_s2 <- list(auc = numeric(5), f1 = numeric(5), brier = numeric(5))

confirmed_interactions_per_fold <- list()

for (f in 1:5) {
    cat(sprintf("\n=== Processing Fold %d/5 ===\n", f))
    train_idx <- which(folds != f)
    test_idx <- which(folds == f)
    
    y_train <- y[train_idx]
    y_test <- y[test_idx]
    
    # 1. Option B: Covariate-only joint imputation (strictly excluding Y)
    cat("  Running missForest (covariates only)...\n")
    set.seed(42 + f)
    imp_res <- missForest(X_miss, ntree = 20)
    X_imp_all <- imp_res$ximp
    
    Z_all <- as.data.frame(lapply(X_miss, function(x) as.numeric(is.na(x))))
    names(Z_all) <- paste0(names(X_imp_all), "_Z")
    
    # Split imputed covariates and Z into train and test
    X_imp_train <- X_imp_all[train_idx, , drop = FALSE]
    Z_train <- Z_all[train_idx, , drop = FALSE]
    
    X_imp_test <- X_imp_all[test_idx, , drop = FALSE]
    Z_test <- Z_all[test_idx, , drop = FALSE]
    
    # Construct base datasets (Main + Indicators)
    X_std_train <- cbind(X_imp_train, Z_train)
    X_std_test <- cbind(X_imp_test, Z_test)
    
    # 2. Run ITIP screening strictly on training data
    cat("  Running ITIP on training fold...\n")
    train_data_for_itip <- data[train_idx, ] # Keep outcome
    set.seed(42 + f)
    # Pass X_imp_train to bypass internal imputation while retaining missing data pattern from train_data_for_itip
    res <- itip(train_data_for_itip, outcome = "mortality", threshold_method = "adaptive", verbose = FALSE, X_imputed = X_imp_train)
    
    # 3. Build Full Interactions dataset
    X_int_train <- X_imp_train
    X_int_test <- X_imp_test
    
    for (col in names(X_imp_train)) {
        z_col <- paste0(col, "_Z")
        if (sum(Z_train[[z_col]]) > 0) {
            # Generate interaction for both train and test (based on Z presence in train)
            X_int_train[[paste0(col, "_int")]] <- X_imp_train[[col]] * Z_train[[z_col]]
            X_int_test[[paste0(col, "_int")]] <- X_imp_test[[col]] * Z_test[[z_col]]
        }
    }
    interaction_cols <- names(X_int_train)[grepl("_int$", names(X_int_train))]
    
    X_full_train <- cbind(X_std_train, X_int_train[, interaction_cols])
    X_full_test <- cbind(X_std_test, X_int_test[, interaction_cols])
    
    # 4. Build Stage 1 Surviving dataset
    stage1_vars <- gsub("^I_", "", res$stage1_kept)
    stage1_ints <- paste0(stage1_vars, "_int")
    # Only keep interactions that were successfully generated
    stage1_ints <- intersect(stage1_ints, interaction_cols)
    
    if (length(stage1_ints) > 0) {
        X_s1_train <- cbind(X_std_train, X_int_train[, stage1_ints, drop=FALSE])
        X_s1_test <- cbind(X_std_test, X_int_test[, stage1_ints, drop=FALSE])
    } else {
        X_s1_train <- X_std_train
        X_s1_test <- X_std_test
    }
    
    # 5. Build Stage 2 Confirmed dataset
    stage2_vars <- gsub("^I_", "", res$pruned_features$kept_interactions)
    stage2_ints <- paste0(stage2_vars, "_int")
    stage2_ints <- intersect(stage2_ints, interaction_cols)
    
    confirmed_interactions_per_fold[[f]] <- stage2_vars
    cat(sprintf("  Confirmed interactions in fold %d: %s\n", f, paste(stage2_vars, collapse=", ")))
    
    if (length(stage2_ints) > 0) {
        X_s2_train <- cbind(X_std_train, X_int_train[, stage2_ints, drop=FALSE])
        X_s2_test <- cbind(X_std_test, X_int_test[, stage2_ints, drop=FALSE])
    } else {
        X_s2_train <- X_std_train
        X_s2_test <- X_std_test
    }
    
    # 6. Evaluate all models
    evaluate_model <- function(X_tr, y_tr, X_te, y_te) {
        set.seed(42 + f)
        fit <- glmnet(as.matrix(X_tr), y_tr, family = "binomial", alpha = 1)
        preds <- predict(fit, as.matrix(X_te), s = min(fit$lambda), type = "response")
        preds <- as.numeric(preds)
        
        auc_val <- auc(y_te, preds, quiet = TRUE)
        y_pred <- as.numeric(preds > 0.5)
        tp <- sum(y_pred == 1 & y_te == 1)
        fp <- sum(y_pred == 1 & y_te == 0)
        fn <- sum(y_pred == 0 & y_te == 1)
        precision <- if (tp + fp > 0) tp / (tp + fp) else 0
        recall <- if (tp + fn > 0) tp / (tp + fn) else 0
        f1_val <- if (precision + recall > 0) 2 * (precision * recall) / (precision + recall) else 0
        brier_val <- mean((preds - y_te)^2)
        
        return(list(auc = auc_val, f1 = f1_val, brier = brier_val))
    }
    
    res_std <- evaluate_model(X_std_train, y_train, X_std_test, y_test)
    metrics_std$auc[f] <- res_std$auc
    metrics_std$f1[f] <- res_std$f1
    metrics_std$brier[f] <- res_std$brier
    
    res_full <- evaluate_model(X_full_train, y_train, X_full_test, y_test)
    metrics_full$auc[f] <- res_full$auc
    metrics_full$f1[f] <- res_full$f1
    metrics_full$brier[f] <- res_full$brier
    
    res_s1 <- evaluate_model(X_s1_train, y_train, X_s1_test, y_test)
    metrics_s1$auc[f] <- res_s1$auc
    metrics_s1$f1[f] <- res_s1$f1
    metrics_s1$brier[f] <- res_s1$brier
    
    res_s2 <- evaluate_model(X_s2_train, y_train, X_s2_test, y_test)
    metrics_s2$auc[f] <- res_s2$auc
    metrics_s2$f1[f] <- res_s2$f1
    metrics_s2$brier[f] <- res_s2$brier
}

cat("\n=== Final Results (5-Fold CV) ===\n")
print_metrics <- function(name, mets) {
    mean_auc <- mean(mets$auc)
    se_auc <- sd(mets$auc) / sqrt(5)
    mean_f1 <- mean(mets$f1)
    se_f1 <- sd(mets$f1) / sqrt(5)
    mean_brier <- mean(mets$brier)
    se_brier <- sd(mets$brier) / sqrt(5)
    
    cat(sprintf("%-20s - AUC: %.4f ± %.4f | F1: %.4f ± %.4f | Brier: %.4f ± %.4f\n",
                name, mean_auc, se_auc, mean_f1, se_f1, mean_brier, se_brier))
}

print_metrics("Standard (Main+Z)", metrics_std)
print_metrics("Full Interactions", metrics_full)
print_metrics("Stage 1 Surviving", metrics_s1)
print_metrics("Stage 2 Confirmed", metrics_s2)

cat("\nInteractions confirmed per fold:\n")
for (f in 1:5) {
    if (length(confirmed_interactions_per_fold[[f]]) == 0) {
        cat(sprintf("Fold %d: None\n", f))
    } else {
        cat(sprintf("Fold %d: %s\n", f, paste(confirmed_interactions_per_fold[[f]], collapse = ", ")))
    }
}

cat("\nPaired t-test results (AUC):\n")
print(t.test(metrics_full$auc, metrics_std$auc, paired=TRUE))

cat("\nPaired t-test results (F1):\n")
print(t.test(metrics_full$f1, metrics_std$f1, paired=TRUE))

cat("\nPaired t-test results (Brier):\n")
print(t.test(metrics_full$brier, metrics_std$brier, paired=TRUE))

