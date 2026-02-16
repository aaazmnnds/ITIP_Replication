# 02_threshold_selection.R
# Methods for selecting optimal pruning threshold ε in ITIP
# Author: Azman Nads
# Supervised by: Daniel Andrade

library(ggplot2)
library(caret) # For cross-validation

if (file.exists("01_ig_calculation.R")) {
  source("01_ig_calculation.R") # Load IG calculation functions
} else {
  source("scripts/01_ig_calculation.R")
}


#' Fixed Threshold Approach
#'
#' Use a predefined threshold value (e.g., ε = 0.01)
#'
#' @param epsilon Fixed threshold value (default: 0.01 = 1% information gain)
#' @return Threshold value
#' @examples
#' epsilon <- select_threshold_fixed(0.01)
select_threshold_fixed <- function(epsilon = 0.01) {
  cat(sprintf("Using fixed threshold: ε = %.3f bits\n", epsilon))
  return(epsilon)
}


#' Adaptive Threshold Based on IG Distribution
#'
#' Set threshold as a fraction of the median IG value
#'
#' @param ig_values Vector of IG values for all interactions
#' @param fraction Fraction of median to use (default: 0.5)
#' @return Adaptive threshold value
#' @examples
#' ig_values <- c(0.45, 0.32, 0.15, 0.08, 0.02, 0.01)
#' epsilon <- select_threshold_adaptive(ig_values)
select_threshold_adaptive <- function(ig_values, fraction = 0.5) {
  # Calculate median IG
  median_ig <- median(ig_values)

  # Set threshold as fraction of median
  epsilon <- median_ig * fraction

  cat(sprintf(
    "Adaptive threshold: ε = %.3f bits (median = %.3f, fraction = %.2f)\n",
    epsilon, median_ig, fraction
  ))

  return(epsilon)
}


#' Percentile-Based Threshold
#'
#' Set threshold at a specific percentile of IG distribution
#'
#' @param ig_values Vector of IG values
#' @param percentile Percentile to use (default: 25th percentile)
#' @return Threshold value
#' @examples
#' ig_values <- c(0.45, 0.32, 0.15, 0.08, 0.02, 0.01)
#' epsilon <- select_threshold_percentile(ig_values, percentile = 25)
select_threshold_percentile <- function(ig_values, percentile = 25) {
  epsilon <- quantile(ig_values, probs = percentile / 100)

  cat(sprintf(
    "Percentile-based threshold: ε = %.3f bits (%dth percentile)\n",
    epsilon, percentile
  ))

  return(epsilon)
}


#' Cross-Validation Threshold Selection
#'
#' Select threshold that maximizes held-out AUC
#'
#' @param data Data frame with outcome and features
#' @param outcome Name of outcome column
#' @param missing_indicators Vector of missing indicator column names
#' @param imputed_vars Vector of imputed variable column names
#' @param epsilon_range Vector of threshold values to test
#' @param n_folds Number of CV folds (default: 5)
#' @return List with optimal threshold and CV results
#' @examples
#' # See example at end of file
select_threshold_cv <- function(data, outcome, missing_indicators, imputed_vars,
                                epsilon_range = seq(0.001, 0.1, by = 0.01),
                                n_folds = 5) {
  cat(sprintf(
    "Running %d-fold cross-validation for %d threshold values...\n",
    n_folds, length(epsilon_range)
  ))

  # Create folds
  set.seed(42)
  folds <- createFolds(data[[outcome]], k = n_folds, list = TRUE)

  # Initialize results
  cv_results <- data.frame(
    epsilon = epsilon_range,
    mean_auc = numeric(length(epsilon_range)),
    sd_auc = numeric(length(epsilon_range))
  )

  # Loop through each threshold
  for (i in seq_along(epsilon_range)) {
    epsilon <- epsilon_range[i]
    auc_folds <- numeric(n_folds)

    # Loop through each fold
    for (fold_idx in seq_along(folds)) {
      # Split data
      test_indices <- folds[[fold_idx]]
      train_data <- data[-test_indices, ]
      test_data <- data[test_indices, ]

      # Calculate IG on training data
      ig_results <- calculate_all_ig(train_data, outcome, missing_indicators, imputed_vars)

      # Prune interactions based on threshold
      keep_vars <- ig_results$variable[ig_results$IG >= epsilon]

      # Build model with pruned features
      # (Simplified: use logistic regression)
      formula_str <- paste(outcome, "~", paste(c(missing_indicators, keep_vars), collapse = " + "))
      model <- glm(as.formula(formula_str), data = train_data, family = binomial())

      # Predict on test data
      predictions <- predict(model, newdata = test_data, type = "response")

      # Calculate AUC
      auc <- pROC::auc(pROC::roc(test_data[[outcome]], predictions, quiet = TRUE))
      auc_folds[fold_idx] <- as.numeric(auc)
    }

    # Store results
    cv_results$mean_auc[i] <- mean(auc_folds)
    cv_results$sd_auc[i] <- sd(auc_folds)
  }

  # Find optimal threshold
  optimal_idx <- which.max(cv_results$mean_auc)
  optimal_epsilon <- cv_results$epsilon[optimal_idx]
  optimal_auc <- cv_results$mean_auc[optimal_idx]

  cat(sprintf(
    "Optimal threshold: ε = %.3f (AUC = %.3f ± %.3f)\n",
    optimal_epsilon, optimal_auc, cv_results$sd_auc[optimal_idx]
  ))

  return(list(
    optimal_epsilon = optimal_epsilon,
    optimal_auc = optimal_auc,
    cv_results = cv_results
  ))
}


#' Visualize Threshold Sensitivity
#'
#' Plot how model performance changes with different threshold values
#'
#' @param cv_results Data frame from select_threshold_cv
#' @return ggplot object
plot_threshold_sensitivity <- function(cv_results) {
  p <- ggplot(cv_results, aes(x = epsilon, y = mean_auc)) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_ribbon(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
      alpha = 0.2, fill = "steelblue"
    ) +
    geom_point(color = "steelblue", size = 2) +
    geom_vline(
      xintercept = cv_results$epsilon[which.max(cv_results$mean_auc)],
      linetype = "dashed", color = "red"
    ) +
    labs(
      title = "Threshold Sensitivity Analysis",
      subtitle = "How pruning threshold affects model performance",
      x = "Pruning Threshold (ε) [bits]",
      y = "Cross-Validated AUC"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12)
    )

  return(p)
}


#' Visualize IG Distribution
#'
#' Plot histogram of IG values with threshold overlay
#'
#' @param ig_values Vector of IG values
#' @param epsilon Threshold value
#' @param variable_names Optional vector of variable names
#' @return ggplot object
plot_ig_distribution <- function(ig_values, epsilon, variable_names = NULL) {
  # Create data frame
  df <- data.frame(IG = ig_values)
  if (!is.null(variable_names)) {
    df$variable <- variable_names
  }

  # Create plot
  p <- ggplot(df, aes(x = IG)) +
    geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = epsilon, linetype = "dashed", color = "red", size = 1) +
    annotate("text",
      x = epsilon, y = Inf, label = sprintf("ε = %.3f", epsilon),
      vjust = 1.5, hjust = -0.1, color = "red", size = 4
    ) +
    labs(
      title = "Information Gain Distribution",
      subtitle = sprintf(
        "%d interactions, %d pruned (IG < %.3f)",
        length(ig_values),
        sum(ig_values < epsilon),
        epsilon
      ),
      x = "Conditional Information Gain (bits)",
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12)
    )

  return(p)
}


#' Compare Threshold Selection Methods
#'
#' Run all threshold selection methods and compare results
#'
#' @param ig_values Vector of IG values
#' @param data Optional data frame for CV method
#' @param outcome Optional outcome column name for CV
#' @param missing_indicators Optional for CV
#' @param imputed_vars Optional for CV
#' @return Data frame comparing methods
compare_threshold_methods <- function(ig_values, data = NULL, outcome = NULL,
                                      missing_indicators = NULL, imputed_vars = NULL) {
  results <- data.frame(
    method = character(),
    epsilon = numeric(),
    n_pruned = numeric(),
    pct_pruned = numeric(),
    stringsAsFactors = FALSE
  )

  # Fixed threshold
  eps_fixed <- select_threshold_fixed(0.01)
  results <- rbind(results, data.frame(
    method = "Fixed (0.01)",
    epsilon = eps_fixed,
    n_pruned = sum(ig_values < eps_fixed),
    pct_pruned = 100 * mean(ig_values < eps_fixed)
  ))

  # Adaptive threshold
  eps_adaptive <- select_threshold_adaptive(ig_values)
  results <- rbind(results, data.frame(
    method = "Adaptive (median/2)",
    epsilon = eps_adaptive,
    n_pruned = sum(ig_values < eps_adaptive),
    pct_pruned = 100 * mean(ig_values < eps_adaptive)
  ))

  # Percentile-based
  eps_percentile <- select_threshold_percentile(ig_values, percentile = 25)
  results <- rbind(results, data.frame(
    method = "Percentile (25th)",
    epsilon = eps_percentile,
    n_pruned = sum(ig_values < eps_percentile),
    pct_pruned = 100 * mean(ig_values < eps_percentile)
  ))

  # Cross-validation (if data provided)
  if (!is.null(data) && !is.null(outcome)) {
    cv_result <- select_threshold_cv(data, outcome, missing_indicators, imputed_vars)
    results <- rbind(results, data.frame(
      method = "Cross-Validation",
      epsilon = cv_result$optimal_epsilon,
      n_pruned = sum(ig_values < cv_result$optimal_epsilon),
      pct_pruned = 100 * mean(ig_values < cv_result$optimal_epsilon)
    ))
  }

  cat("\n--- Threshold Comparison ---\n")
  print(results)

  return(results)
}


#' Statistical Significance Test for IG
#'
#' Test if IG is significantly greater than 0 using permutation test
#'
#' @param Y Outcome vector
#' @param Z Missing indicator
#' @param I Interaction term
#' @param n_permutations Number of permutations (default: 1000)
#' @param alpha Significance level (default: 0.05)
#' @return List with p-value and significance decision
test_ig_significance <- function(Y, Z, I, n_permutations = 1000, alpha = 0.05) {
  # Calculate observed IG
  IG_observed <- calculate_conditional_ig_discrete(Y, Z, I)

  # Permutation test: shuffle Y to break association
  IG_permuted <- numeric(n_permutations)

  for (i in 1:n_permutations) {
    Y_perm <- sample(Y)
    IG_permuted[i] <- calculate_conditional_ig_discrete(Y_perm, Z, I)
  }

  # Calculate p-value
  p_value <- mean(IG_permuted >= IG_observed)

  # Significance decision
  is_significant <- p_value < alpha

  cat(sprintf(
    "IG = %.3f, p-value = %.3f, significant = %s (α = %.2f)\n",
    IG_observed, p_value, is_significant, alpha
  ))

  return(list(
    IG = IG_observed,
    p_value = p_value,
    is_significant = is_significant,
    permuted_distribution = IG_permuted
  ))
}


# Example usage (for testing)
if (FALSE) {
  # Simulate IG values
  set.seed(42)
  ig_values <- c(
    runif(5, 0.3, 0.5), # High IG interactions
    runif(10, 0.05, 0.2), # Medium IG interactions
    runif(15, 0.001, 0.05) # Low IG interactions
  )

  # Compare threshold methods
  cat("--- Threshold Selection Comparison ---\n")
  compare_threshold_methods(ig_values)

  # Visualize IG distribution
  epsilon <- select_threshold_fixed(0.01)
  plot_ig_distribution(ig_values, epsilon)

  # Test significance
  Y <- sample(c(0, 1), 100, replace = TRUE)
  Z <- sample(c(0, 1), 100, replace = TRUE, prob = c(0.7, 0.3))
  I <- rnorm(100) * Z
  test_ig_significance(Y, Z, I)
}
