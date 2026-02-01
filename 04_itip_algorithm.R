# 04_itip_algorithm.R
# Main ITIP (Information-Theoretic Interaction Pruning) Algorithm
# Author: Azman Nads
# Supervised by: Daniel Andrade

library(missForest) # For imputation
library(BAS) # For Bayesian Adaptive Sampling
library(dplyr) # For data manipulation

source("01_ig_calculation.R")
source("02_threshold_selection.R")
source("03_bas_utils.R")

#' ITIP: Information-Theoretic Interaction Pruning
#'
#' Main algorithm that implements the complete ITIP workflow:
#' 1. Identify missing data patterns
#' 2. Create missing indicators (Z) - BEFORE imputation
#' 3. Impute missing data using missForest
#' 4. Generate interaction terms (X_imp × Z)
#' 5. Calculate Conditional Information Gain for each interaction
#' 6. Prune interactions where IG < epsilon
#' 7. Return pruned feature set
#'
#' @param data Data frame with missing values
#' @param outcome Name of outcome column (binary)
#' @param epsilon Pruning threshold (default: 0.01)
#' @param threshold_method Method for selecting epsilon: "fixed", "adaptive", "cv"
#' @param verbose Print progress messages (default: TRUE)
#' @return List containing:
#'   - data_imputed: Imputed dataset
#'   - missing_indicators: Binary indicators for missingness
#'   - interactions: Generated interaction terms
#'   - ig_results: IG values for all interactions
#'   - pruned_features: Features to keep after pruning
#'   - bas_results: BAS model results (PIP, best model)
#'   - pruning_stats: Summary statistics
#' @examples
#' # See examples at end of file
itip <- function(data, outcome, epsilon = 0.01,
                 threshold_method = "fixed", verbose = TRUE) {
  if (verbose) cat("=== ITIP Algorithm ===\n")

  # Step 1: Identify variables with missing data
  if (verbose) cat("Step 1: Identifying missing data patterns...\n")

  # Separate outcome from features
  Y <- data[[outcome]]
  X <- data[, setdiff(names(data), outcome), drop = FALSE]

  # Identify columns with missing values
  missing_cols <- names(X)[colSums(is.na(X)) > 0]

  if (length(missing_cols) == 0) {
    warning("No missing data found. ITIP is not needed.")
    return(list(
      data_imputed = data,
      missing_indicators = NULL,
      interactions = NULL,
      ig_results = NULL,
      pruned_features = list(
        imputed_vars = names(X),
        missing_indicators = character(0),
        kept_interactions = character(0)
      ),
      pruning_stats = list(n_pruned = 0, pct_pruned = 0)
    ))
  }

  if (verbose) {
    cat(sprintf("  Found %d variables with missing data\n", length(missing_cols)))
    cat(sprintf("  Missing percentages:\n"))
    for (col in missing_cols) {
      pct_missing <- 100 * mean(is.na(X[[col]]))
      cat(sprintf("    %s: %.1f%%\n", col, pct_missing))
    }
  }

  # Step 2: Create missing indicators BEFORE imputation
  if (verbose) cat("\nStep 2: Creating missing indicators (before imputation)...\n")

  missing_indicators <- data.frame(matrix(0, nrow = nrow(X), ncol = length(missing_cols)))
  names(missing_indicators) <- paste0("Z_", missing_cols)

  for (i in seq_along(missing_cols)) {
    col <- missing_cols[i]
    missing_indicators[[paste0("Z_", col)]] <- as.integer(is.na(X[[col]]))
  }

  if (verbose) {
    cat(sprintf("  Created %d missing indicators\n", ncol(missing_indicators)))
  }

  # Step 3: Impute missing data using missForest
  if (verbose) cat("\nStep 3: Imputing missing data with missForest...\n")

  # missForest requires all numeric or all factor
  # Convert to appropriate types
  X_for_imputation <- X

  # Run missForest
  set.seed(42) # For reproducibility
  imputation_result <- missForest(X_for_imputation, verbose = verbose)

  X_imputed <- imputation_result$ximp

  if (verbose) {
    cat(sprintf("  Imputation error (OOB): %.4f\n", imputation_result$OOBerror))
  }

  # Step 4: Generate interaction terms
  if (verbose) cat("\nStep 4: Generating interaction terms (X_imp × Z)...\n")

  interactions <- data.frame(matrix(0, nrow = nrow(X), ncol = length(missing_cols)))
  names(interactions) <- paste0("I_", missing_cols)

  for (i in seq_along(missing_cols)) {
    col <- missing_cols[i]
    Z_col <- paste0("Z_", col)
    I_col <- paste0("I_", col)

    # Interaction: X_imputed × Z
    interactions[[I_col]] <- X_imputed[[col]] * missing_indicators[[Z_col]]
  }

  if (verbose) {
    cat(sprintf("  Generated %d interaction terms\n", ncol(interactions)))
  }

  # Step 5: Calculate Conditional Information Gain
  if (verbose) cat("\nStep 5: Calculating Conditional Information Gain...\n")

  ig_results <- data.frame(
    variable = character(),
    IG = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(missing_cols)) {
    col <- missing_cols[i]
    Z <- missing_indicators[[paste0("Z_", col)]]
    I <- interactions[[paste0("I_", col)]]

    # Calculate IG(I | Z)
    IG <- calculate_conditional_ig_discrete(Y, Z, I)

    ig_results <- rbind(ig_results, data.frame(
      variable = col,
      IG = IG,
      stringsAsFactors = FALSE
    ))

    if (verbose) {
      cat(sprintf("  %s: IG = %.4f bits\n", col, IG))
    }
  }

  # Sort by IG (descending)
  ig_results <- ig_results[order(-ig_results$IG), ]

  # Step 6: Select pruning threshold
  if (verbose) cat("\nStep 6: Selecting pruning threshold...\n")

  if (threshold_method == "fixed") {
    # epsilon already provided
    if (verbose) cat(sprintf("  Using fixed threshold: ε = %.4f bits\n", epsilon))
  } else if (threshold_method == "adaptive") {
    epsilon <- select_threshold_adaptive(ig_results$IG)
  } else if (threshold_method == "percentile") {
    epsilon <- select_threshold_percentile(ig_results$IG, percentile = 25)
  } else {
    warning("Unknown threshold method. Using fixed threshold.")
  }

  # Step 7: Prune interactions
  if (verbose) cat("\nStep 7: Pruning low-IG interactions...\n")

  # Identify interactions to keep
  keep_interactions <- ig_results$variable[ig_results$IG >= epsilon]
  prune_interactions <- ig_results$variable[ig_results$IG < epsilon]

  n_pruned <- length(prune_interactions)
  pct_pruned <- 100 * n_pruned / length(missing_cols)

  if (verbose) {
    cat(sprintf(
      "  Pruned %d / %d interactions (%.1f%%)\n",
      n_pruned, length(missing_cols), pct_pruned
    ))
    if (n_pruned > 0) {
      cat("  Pruned interactions:\n")
      for (var in prune_interactions) {
        ig_val <- ig_results$IG[ig_results$variable == var]
        cat(sprintf("    %s (IG = %.4f < %.4f)\n", var, ig_val, epsilon))
      }
    }
    if (length(keep_interactions) > 0) {
      cat("  Kept interactions:\n")
      for (var in keep_interactions) {
        ig_val <- ig_results$IG[ig_results$variable == var]
        cat(sprintf("    %s (IG = %.4f >= %.4f)\n", var, ig_val, epsilon))
      }
    }
  }

  # Step 8: Construct final feature set
  if (verbose) cat("\nStep 8: Constructing final feature set...\n")

  # Always include:
  # 1. All imputed features (X_imputed)
  # 2. All missing indicators (Z)
  # 3. Only high-IG interactions (I where IG >= epsilon)

  pruned_features <- list(
    imputed_vars = names(X_imputed),
    missing_indicators = names(missing_indicators),
    kept_interactions = paste0("I_", keep_interactions)
  )

  if (verbose) {
    cat(sprintf("  Final feature set:\n"))
    cat(sprintf("    - %d imputed variables\n", length(pruned_features$imputed_vars)))
    cat(sprintf("    - %d missing indicators\n", length(pruned_features$missing_indicators)))
    cat(sprintf(
      "    - %d interactions (pruned %d)\n",
      length(pruned_features$kept_interactions), n_pruned
    ))
    cat(sprintf(
      "  Total features: %d\n",
      length(pruned_features$imputed_vars) +
        length(pruned_features$missing_indicators) +
        length(pruned_features$kept_interactions)
    ))
  }

  # Combine all data
  data_imputed <- cbind(
    X_imputed,
    missing_indicators,
    interactions[, pruned_features$kept_interactions, drop = FALSE]
  )
  data_imputed[[outcome]] <- Y

  # Step 9: Run BAS on final feature set
  if (verbose) cat("\nStep 9: Running Bayesian Adaptive Sampling (BAS)...\n")

  # Construct formula for BAS
  bas_formula_str <- paste(outcome, "~ .")
  bas_formula <- as.formula(bas_formula_str)

  # Run BAS
  bas_results <- tryCatch(
    {
      run_bas_selection(
        formula = bas_formula,
        data = data_imputed,
        method = "MCMC",
        verbose = verbose
      )
    },
    error = function(e) {
      if (verbose) cat(sprintf("  Warning: BAS failed (%s). Returning NULL.\n", e$message))
      return(NULL)
    }
  )

  if (verbose) cat("\n=== ITIP Complete ===\n")

  # Return results
  return(list(
    data_imputed = data_imputed,
    missing_indicators = missing_indicators,
    interactions = interactions,
    ig_results = ig_results,
    pruned_features = pruned_features,
    bas_results = bas_results,
    epsilon = epsilon,
    epsilon = epsilon,
    pruning_stats = list(
      n_total = length(missing_cols),
      n_pruned = n_pruned,
      pct_pruned = pct_pruned,
      n_kept = length(keep_interactions)
    )
  ))
}


#' Build prediction model with ITIP-pruned features
#'
#' @param itip_result Result from itip() function
#' @param outcome Name of outcome column
#' @param method Modeling method: "logistic", "rf", "xgboost"
#' @return Fitted model
build_itip_model <- function(itip_result, outcome, method = "logistic") {
  data <- itip_result$data_imputed

  # Create formula
  feature_names <- c(
    itip_result$pruned_features$imputed_vars,
    itip_result$pruned_features$missing_indicators,
    itip_result$pruned_features$kept_interactions
  )

  formula_str <- paste(outcome, "~", paste(feature_names, collapse = " + "))
  formula_obj <- as.formula(formula_str)

  if (method == "logistic") {
    model <- glm(formula_obj, data = data, family = binomial())
  } else if (method == "rf") {
    library(randomForest)
    model <- randomForest(formula_obj, data = data)
  } else {
    stop("Unknown method. Use 'logistic' or 'rf'")
  }

  return(model)
}


#' Compare ITIP vs. Baseline (all interactions)
#'
#' @param data Data frame with missing values
#' @param outcome Name of outcome column
#' @param epsilon Pruning threshold
#' @param n_folds Number of CV folds (default: 5)
#' @return Comparison results
compare_itip_baseline <- function(data, outcome, epsilon = 0.01, n_folds = 5) {
  cat("=== Comparing ITIP vs. Baseline ===\n\n")

  library(caret)
  library(pROC)

  # Create folds
  set.seed(42)
  folds <- createFolds(data[[outcome]], k = n_folds, list = TRUE)

  results <- data.frame(
    fold = integer(),
    method = character(),
    auc = numeric(),
    n_features = integer(),
    stringsAsFactors = FALSE
  )

  for (fold_idx in seq_along(folds)) {
    cat(sprintf("Fold %d/%d\n", fold_idx, n_folds))

    # Split data
    test_indices <- folds[[fold_idx]]
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]

    # ITIP model
    itip_result <- itip(train_data, outcome, epsilon, verbose = FALSE)
    itip_model <- build_itip_model(itip_result, outcome, method = "logistic")

    # Predict on test data
    test_data_itip <- test_data
    # Need to apply same imputation and feature engineering...
    # (Simplified for now - in practice, need to handle this carefully)

    # For now, just record feature counts
    n_features_itip <- length(itip_result$pruned_features$imputed_vars) +
      length(itip_result$pruned_features$missing_indicators) +
      length(itip_result$pruned_features$kept_interactions)

    cat(sprintf(
      "  ITIP: %d features (pruned %d interactions)\n",
      n_features_itip, itip_result$pruning_stats$n_pruned
    ))
  }

  cat("\n=== Comparison Complete ===\n")

  return(results)
}


# Example usage (for testing)
if (FALSE) {
  # Create toy dataset
  set.seed(42)
  n <- 200

  data <- data.frame(
    X1 = rnorm(n),
    X2 = rnorm(n),
    X3 = rnorm(n),
    Y = sample(c(0, 1), n, replace = TRUE)
  )

  # Introduce missing data
  data$X1[sample(1:n, 30)] <- NA
  data$X2[sample(1:n, 40)] <- NA

  # Run ITIP
  result <- itip(data, outcome = "Y", epsilon = 0.01, verbose = TRUE)

  # Inspect results
  print(result$ig_results)
  print(result$pruning_stats)
}
