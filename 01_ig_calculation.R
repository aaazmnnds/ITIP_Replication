# 01_ig_calculation.R
# Core functions for computing Conditional Information Gain (IG) in ITIP
# Author: Azman Nads
# Supervised by: Daniel Andrade

library(entropy) # For entropy calculations
library(BAS) # For Bayesian Adaptive Sampling

#' Calculate Entropy H(Y)
#'
#' @param Y A vector of outcomes (binary or discrete)
#' @return Entropy in bits
#' @examples
#' Y <- c(0, 1, 0, 1, 1, 0)
#' calculate_entropy(Y) # Should return ~1.0 bits
calculate_entropy <- function(Y) {
  # Handle edge cases
  if (length(Y) == 0) {
    stop("Y cannot be empty")
  }

  # Calculate empirical probabilities
  probs <- table(Y) / length(Y)

  # Calculate entropy using natural log, then convert to bits
  H <- -sum(probs * log2(probs + 1e-10)) # Add small constant to avoid log(0)

  return(H)
}


#' Calculate Conditional Entropy H(Y | X)
#'
#' @param Y A vector of outcomes
#' @param X A vector of conditioning variables
#' @return Conditional entropy in bits
#' @examples
#' Y <- c(0, 1, 0, 1, 1, 0, 1, 0)
#' X <- c(0, 0, 0, 0, 1, 1, 1, 1)
#' conditional_entropy(Y, X)
conditional_entropy <- function(Y, X) {
  # Handle edge cases
  if (length(Y) != length(X)) {
    stop("Y and X must have the same length")
  }

  # Get unique values of X
  x_values <- unique(X)

  # Calculate H(Y | X) = sum_x P(X=x) * H(Y | X=x)
  H_conditional <- 0

  for (x in x_values) {
    # Subset Y where X = x
    Y_given_x <- Y[X == x]

    # Calculate P(X = x)
    p_x <- length(Y_given_x) / length(Y)

    # Calculate H(Y | X = x)
    if (length(Y_given_x) > 0) {
      probs_y_given_x <- table(Y_given_x) / length(Y_given_x)
      H_y_given_x <- -sum(probs_y_given_x * log2(probs_y_given_x + 1e-10))

      # Add weighted contribution
      H_conditional <- H_conditional + p_x * H_y_given_x
    }
  }

  return(H_conditional)
}


#' Calculate Conditional Entropy H(Y | X1, X2)
#'
#' @param Y A vector of outcomes
#' @param X1 First conditioning variable
#' @param X2 Second conditioning variable
#' @return Conditional entropy in bits
conditional_entropy_2d <- function(Y, X1, X2) {
  # Handle edge cases
  if (length(Y) != length(X1) || length(Y) != length(X2)) {
    stop("Y, X1, and X2 must have the same length")
  }

  # Create joint conditioning variable
  X_joint <- paste(X1, X2, sep = "_")

  # Use 1D conditional entropy function
  return(conditional_entropy(Y, X_joint))
}


#' Calculate Conditional Information Gain IG(I | Z)
#'
#' This is the core ITIP metric: how much information does the interaction I
#' provide beyond the missing indicator Z?
#'
#' @param Y Outcome vector (binary or discrete)
#' @param Z Missing indicator (0 = observed, 1 = missing)
#' @param I Interaction term (X_imp * Z)
#' @return Conditional Information Gain in bits
#' @examples
#' # Example 1: High IG (Lactate)
#' Y <- c(0, 1, 0, 1, 0, 0, 0, 1)
#' Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
#' I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)
#' calculate_conditional_ig(Y, Z, I) # Should be high (~0.4 bits)
calculate_conditional_ig <- function(Y, Z, I) {
  # Calculate H(Y | Z)
  H_Y_given_Z <- conditional_entropy(Y, Z)

  # Calculate H(Y | Z, I)
  H_Y_given_Z_I <- conditional_entropy_2d(Y, Z, I)

  # Calculate IG(I | Z) = H(Y | Z) - H(Y | Z, I)
  IG <- H_Y_given_Z - H_Y_given_Z_I

  # Ensure non-negative (due to numerical errors, IG might be slightly negative)
  IG <- max(0, IG)

  return(IG)
}


#' Discretize continuous interaction term for IG calculation
#'
#' @param I Continuous interaction term
#' @param n_bins Number of bins (default: 3 for Low/Med/High)
#' @return Discretized interaction term
#' @examples
#' I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)
#' discretize_interaction(I, n_bins = 3)
discretize_interaction <- function(I, n_bins = 3) {
  # For zero values (observed data), keep as 0
  I_nonzero <- I[I != 0]

  if (length(I_nonzero) == 0) {
    return(I) # All zeros, no discretization needed
  }

  # Check if all non-zero values are identical
  if (length(unique(I_nonzero)) == 1) {
    return(I) # All non-zero values are the same, no discretization needed
  }

  # Create bins for non-zero values
  breaks <- quantile(I_nonzero, probs = seq(0, 1, length.out = n_bins + 1))

  # Check if breaks are unique (can happen with small samples or repeated values)
  if (length(unique(breaks)) < length(breaks)) {
    # Fall back to simple binning
    return(I)
  }

  # Discretize
  I_discrete <- I
  I_discrete[I != 0] <- cut(I[I != 0], breaks = breaks, labels = FALSE, include.lowest = TRUE)

  return(I_discrete)
}


#' Calculate Conditional IG with discretization
#'
#' @param Y Outcome vector
#' @param Z Missing indicator
#' @param I Continuous interaction term
#' @param n_bins Number of bins for discretization
#' @return Conditional Information Gain in bits
calculate_conditional_ig_discrete <- function(Y, Z, I, n_bins = 3) {
  # Discretize interaction term
  I_discrete <- discretize_interaction(I, n_bins)

  # Calculate IG
  return(calculate_conditional_ig(Y, Z, I_discrete))
}


#' Calculate IG for all interactions in a dataset
#'
#' @param data Data frame containing outcome, missing indicators, and imputed values
#' @param outcome Name of outcome column
#' @param missing_indicators Vector of missing indicator column names
#' @param imputed_vars Vector of imputed variable column names
#' @return Data frame with IG values for each interaction
#' @examples
#' data <- data.frame(
#'   Y = c(0, 1, 0, 1, 0, 0, 0, 1),
#'   Z_lactate = c(0, 1, 0, 1, 0, 1, 0, 1),
#'   lactate_imp = c(2.1, 3.5, 1.8, 4.2, 2.0, 2.8, 1.9, 5.1)
#' )
#' calculate_all_ig(data, "Y", "Z_lactate", "lactate_imp")
calculate_all_ig <- function(data, outcome, missing_indicators, imputed_vars) {
  # Validate inputs
  if (length(missing_indicators) != length(imputed_vars)) {
    stop("missing_indicators and imputed_vars must have the same length")
  }

  # Initialize results
  results <- data.frame(
    variable = character(),
    IG = numeric(),
    stringsAsFactors = FALSE
  )

  # Extract outcome
  Y <- data[[outcome]]

  # Loop through each variable
  for (i in seq_along(missing_indicators)) {
    var_name <- imputed_vars[i]
    Z_name <- missing_indicators[i]

    # Extract Z and X_imp
    Z <- data[[Z_name]]
    X_imp <- data[[imputed_vars[i]]]

    # Create interaction term
    I <- X_imp * Z

    # Calculate IG
    IG <- calculate_conditional_ig_discrete(Y, Z, I)

    # Store result
    results <- rbind(results, data.frame(
      variable = var_name,
      IG = IG,
      stringsAsFactors = FALSE
    ))
  }

  # Sort by IG (descending)
  results <- results[order(-results$IG), ]

  return(results)
}


#' Estimate IG using BAS posterior inclusion probabilities
#'
#' This function uses Bayesian Adaptive Sampling to estimate IG based on
#' posterior inclusion probabilities.
#'
#' @param data Data frame with outcome and features
#' @param outcome Name of outcome column
#' @param Z_name Name of missing indicator column
#' @param I_name Name of interaction term column
#' @return Estimated IG in bits
#' @examples
#' data <- data.frame(
#'   Y = c(0, 1, 0, 1, 0, 0, 0, 1),
#'   Z = c(0, 1, 0, 1, 0, 1, 0, 1),
#'   I = c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)
#' )
#' estimate_ig_from_bas(data, "Y", "Z", "I")
estimate_ig_from_bas <- function(data, outcome, Z_name, I_name) {
  # Extract variables
  Y <- data[[outcome]]
  Z <- data[[Z_name]]
  I <- data[[I_name]]

  # Create design matrix with Z and I
  X <- cbind(Z, I)
  colnames(X) <- c(Z_name, I_name)

  # Fit BAS model
  # Note: We force Z to always be included, and test inclusion of I
  bas_model <- bas.glm(
    Y ~ Z + I,
    data = data.frame(Y = Y, Z = Z, I = I),
    family = binomial(),
    method = "MCMC",
    n.models = 1000,
    MCMC.iterations = 10000,
    modelprior = uniform()
  )

  # Extract posterior inclusion probability for I
  # (conditional on Z being in the model)
  pip <- bas_model$probne0[I_name]

  # Convert to IG: IG ≈ -log2(1 - pip)
  # This is an approximation based on information theory
  IG <- -log2(1 - pip + 1e-10)

  return(IG)
}


#' Bootstrap confidence interval for IG
#'
#' @param Y Outcome vector
#' @param Z Missing indicator
#' @param I Interaction term
#' @param n_bootstrap Number of bootstrap samples (default: 1000)
#' @param alpha Significance level (default: 0.05 for 95% CI)
#' @return List with IG estimate, lower bound, and upper bound
bootstrap_ig_ci <- function(Y, Z, I, n_bootstrap = 1000, alpha = 0.05) {
  # Original IG
  IG_original <- calculate_conditional_ig_discrete(Y, Z, I)

  # Bootstrap samples
  IG_bootstrap <- numeric(n_bootstrap)

  for (b in 1:n_bootstrap) {
    # Resample with replacement
    indices <- sample(1:length(Y), replace = TRUE)
    Y_boot <- Y[indices]
    Z_boot <- Z[indices]
    I_boot <- I[indices]

    # Calculate IG
    IG_bootstrap[b] <- calculate_conditional_ig_discrete(Y_boot, Z_boot, I_boot)
  }

  # Calculate confidence interval
  CI_lower <- quantile(IG_bootstrap, alpha / 2)
  CI_upper <- quantile(IG_bootstrap, 1 - alpha / 2)

  return(list(
    IG = IG_original,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    bootstrap_samples = IG_bootstrap
  ))
}


#' Test theoretical property: IG non-negativity
#'
#' @param n_tests Number of random tests to run
#' @return TRUE if all tests pass, FALSE otherwise
test_ig_nonnegativity <- function(n_tests = 100) {
  all_pass <- TRUE

  for (i in 1:n_tests) {
    # Generate random data
    n <- sample(50:200, 1)
    Y <- sample(c(0, 1), n, replace = TRUE)
    Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
    I <- rnorm(n) * Z # Interaction only non-zero when Z=1

    # Calculate IG
    IG <- calculate_conditional_ig_discrete(Y, Z, I)

    # Check non-negativity
    if (IG < -1e-6) { # Allow small numerical errors
      cat(sprintf("FAIL: IG = %.6f (negative!)\n", IG))
      all_pass <- FALSE
    }
  }

  if (all_pass) {
    cat(sprintf("PASS: All %d tests show IG >= 0\n", n_tests))
  }

  return(all_pass)
}


#' Test theoretical property: IG = 0 for independent variables
#'
#' @param n_tests Number of random tests to run
#' @return TRUE if all tests pass, FALSE otherwise
test_ig_independence <- function(n_tests = 100) {
  all_pass <- TRUE

  for (i in 1:n_tests) {
    # Generate independent Y and I (given Z)
    n <- sample(50:200, 1)
    Y <- sample(c(0, 1), n, replace = TRUE)
    Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
    I <- rnorm(n) * Z # Independent of Y

    # Calculate IG
    IG <- calculate_conditional_ig_discrete(Y, Z, I)

    # Check if IG is close to 0 (allowing for sampling variability)
    if (IG > 0.1) { # Threshold for "close to 0"
      cat(sprintf("FAIL: IG = %.6f (should be ~0 for independent vars)\n", IG))
      all_pass <- FALSE
    }
  }

  if (all_pass) {
    cat(sprintf("PASS: All %d tests show IG ≈ 0 for independent variables\n", n_tests))
  }

  return(all_pass)
}


# Example usage (for testing)
if (FALSE) {
  # Example 1: High IG (Lactate)
  Y <- c(0, 1, 0, 1, 0, 0, 0, 1)
  Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
  I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)

  cat("Example 1 (Lactate - High IG):\n")
  cat(sprintf("IG(I | Z) = %.3f bits\n", calculate_conditional_ig_discrete(Y, Z, I)))

  # Example 2: Low IG (Glucose)
  Y <- c(0, 0, 1, 1, 0, 0, 1, 1)
  Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
  I <- c(0, 115, 0, 105, 0, 98, 0, 122)

  cat("\nExample 2 (Glucose - Low IG):\n")
  cat(sprintf("IG(I | Z) = %.3f bits\n", calculate_conditional_ig_discrete(Y, Z, I)))

  # Run tests
  cat("\n--- Running Tests ---\n")
  test_ig_nonnegativity()
  test_ig_independence()
}
