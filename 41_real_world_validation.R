# real_world_validation.R
# Applying ITIP to MIMIC Septic Shock Mortality Prediction

library(BAS)
library(missForest)
library(glmnet)
library(pROC)

# 1. Load Data
cat("Loading MIMIC Septic Shock Data...\n")
data <- read.csv("../data/mimic_septic_shock_final.csv")
# Remove HADM_ID
data <- data[, !(names(data) %in% c("HADM_ID"))]

# Sample for speed during testing if needed (but 2k rows should be fine)
# data <- data[sample(nrow(data), 1000), ]

# Outcome
y <- data$mortality
X_miss <- data[, names(data) != "mortality"]

# 2. Impute with missForest
cat("Imputing missing values with missForest...\n")
# Parallelize if possible, but 2k rows is small
imp_res <- missForest(X_miss, ntree = 20)
X_imp <- imp_res$ximp

# 3. Create Indicators and Interactions
cat("Creating missing indicators and interaction terms...\n")
Z <- as.data.frame(lapply(X_miss, function(x) as.numeric(is.na(x))))
names(Z) <- paste0(names(X_miss), "_Z")

# Interaction terms (only for columns that had missing data AND were kept by missForest)
X_int <- X_imp
kept_cols <- names(X_imp)
Z <- Z[, paste0(kept_cols, "_Z")] # Sync Z with X_imp

for (col in kept_cols) {
    z_col <- paste0(col, "_Z")
    if (sum(Z[[z_col]]) > 0) {
        X_int[[paste0(col, "_int")]] <- X_imp[[col]] * Z[[z_col]]
    }
}

# Add indicators to the feature set
# Proper ITIP model includes: X_imp, Z, and filtered interactions
X_plus_Z <- cbind(X_imp, Z)

# 4. ITIP Pruning
cat("Starting ITIP Pruning (Conditional Information Gain)...\n")
# We loop through each interaction and compute PIP given the indicator
selected_interactions <- c()
interaction_cols <- names(X_int)[grepl("_int$", names(X_int))]

# To speed up, we can use a simpler BAS call per interaction group
# or a full model and check PIPs. JMLR suggests PIP is a good proxy.
# We'll use the scale-aligned proxy: -log(1-gamma)

# For each interaction, we need the model: Y ~ Z_j + Interaction_j
# This isoloates the conditional contribution.
ig_results <- data.frame(Interaction = interaction_cols, IG = 0, PIP = 0)

for (i in seq_along(interaction_cols)) {
    int_col <- interaction_cols[i]
    base_col <- gsub("_int$", "", int_col)
    z_col <- paste0(base_col, "_Z")

    # Model: Y ~ Z + Int
    temp_df <- data.frame(y = y, z = Z[[z_col]], int = X_int[[int_col]])
    # Use BAS to get PIP for 'int'
    bas_fit <- bas.glm(y ~ z + int,
        data = temp_df, family = binomial(),
        include.always = ~z,
        method = "BAS", update = 500, prob.rw = 0.5
    )

    pip <- bas_fit$probne0[3] # Index 3 is 'int'
    ig_results$PIP[i] <- pip
    ig_results$IG[i] <- -log(1 - min(pip, 0.999))

    if (i %% 5 == 0) cat(sprintf("  Processed %d/%d interactions...\n", i, length(interaction_cols)))
}

# Threshold selection (Adaptive: Median/2)
threshold <- median(ig_results$IG) / 2
retained_ints <- ig_results$Interaction[ig_results$IG >= threshold]
cat(sprintf("Retained %d out of %d interactions.\n", length(retained_ints), length(interaction_cols)))

# 5. Final Model Comparison
cat("Evaluating models (5-fold CV)...\n")
set.seed(42)
folds <- sample(rep(1:5, length.out = nrow(X_imp)))

evaluate_set <- function(X_set, y) {
    aucs <- numeric(5)
    f1s <- numeric(5)
    briers <- numeric(5)

    for (f in 1:5) {
        train_idx <- which(folds != f)
        test_idx <- which(folds == f)

        # Fit model
        fit <- glmnet(as.matrix(X_set[train_idx, ]), y[train_idx], family = "binomial", alpha = 1)

        # Select lambda via internal CV or just use a small one
        preds <- predict(fit, as.matrix(X_set[test_idx, ]), s = min(fit$lambda), type = "response")
        preds <- as.numeric(preds)

        # AUC
        aucs[f] <- auc(y[test_idx], preds, quiet = TRUE)

        # F1 (at 0.5 threshold)
        y_pred <- as.numeric(preds > 0.5)
        y_true <- y[test_idx]
        tp <- sum(y_pred == 1 & y_true == 1)
        fp <- sum(y_pred == 1 & y_true == 0)
        fn <- sum(y_pred == 0 & y_true == 1)

        precision <- if (tp + fp > 0) tp / (tp + fp) else 0
        recall <- if (tp + fn > 0) tp / (tp + fn) else 0
        f1s[f] <- if (precision + recall > 0) 2 * (precision * recall) / (precision + recall) else 0

        # Brier Score
        briers[f] <- mean((preds - y_true)^2)
    }
    return(c(AUC = mean(aucs), F1 = mean(f1s), Brier = mean(briers)))
}

# Model 1: Standard (No interactions)
cat("Training Standard Model...\n")
res_std <- evaluate_set(X_plus_Z, y)
cat(sprintf("  Standard - AUC: %.4f, F1: %.4f, Brier: %.4f\n", res_std[1], res_std[2], res_std[3]))

# Model 2: Full Interactions
cat("Training Full Interactions Model...\n")
X_full_int <- cbind(X_plus_Z, X_int[, interaction_cols])
res_full <- evaluate_set(X_full_int, y)
cat(sprintf("  Full Int - AUC: %.4f, F1: %.4f, Brier: %.4f\n", res_full[1], res_full[2], res_full[3]))

# Model 3: ITIP Pruned
cat("Training ITIP Pruned Model...\n")
X_itip <- cbind(X_plus_Z, X_int[, retained_ints])
res_itip <- evaluate_set(X_itip, y)
cat(sprintf("  ITIP     - AUC: %.4f, F1: %.4f, Brier: %.4f\n", res_itip[1], res_itip[2], res_itip[3]))

# Save Results
res_table <- data.frame(
    Model = c("Standard", "Full Interactions", "ITIP"),
    AUC = c(res_std[1], res_full[1], res_itip[1]),
    F1 = c(res_std[2], res_full[2], res_itip[2]),
    Brier = c(res_std[3], res_full[3], res_itip[3])
)
write.csv(res_table, "../results/mimic_validation_results.csv", row.names = FALSE)
write.csv(ig_results, "../results/mimic_ig_pips.csv", row.names = FALSE)
