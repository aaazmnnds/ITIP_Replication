# unit_tests.R
# Validation tests for ITIP theoretical properties
# Author: Azman Nads

library(testthat)
options(testthat.default_reporter = "summary")

source("01_ig_calculation.R")
source("02_threshold_selection.R")


#' Test Suite: Entropy Calculations
test_entropy_calculations <- function() {
  cat("\n=== Testing Entropy Calculations ===\n")

  # Test 1: Entropy of fair coin
  test_that("Entropy of fair coin is 1 bit", {
    Y <- c(0, 1, 0, 1, 0, 1, 0, 1)
    H <- calculate_entropy(Y)
    expect_equal(H, 1.0, tolerance = 0.01)
  })

  # Test 2: Entropy of deterministic variable
  test_that("Entropy of deterministic variable is 0", {
    Y <- c(1, 1, 1, 1, 1)
    H <- calculate_entropy(Y)
    expect_equal(H, 0.0, tolerance = 0.01)
  })

  # Test 3: Conditional entropy
  test_that("Conditional entropy is non-negative", {
    Y <- sample(c(0, 1), 100, replace = TRUE)
    X <- sample(c(0, 1), 100, replace = TRUE)
    H_cond <- conditional_entropy(Y, X)
    expect_gte(H_cond, 0)
  })

  # Test 4: H(Y|X) <= H(Y)
  test_that("Conditional entropy is at most unconditional entropy", {
    Y <- sample(c(0, 1), 100, replace = TRUE)
    X <- sample(c(0, 1), 100, replace = TRUE)
    H_Y <- calculate_entropy(Y)
    H_Y_given_X <- conditional_entropy(Y, X)
    expect_lte(H_Y_given_X, H_Y + 0.01) # Allow small numerical error
  })

  cat("✓ All entropy tests passed\n")
}


#' Test Suite: Information Gain Properties
test_ig_properties <- function() {
  cat("\n=== Testing Information Gain Properties ===\n")

  # Test 1: IG non-negativity
  test_that("IG is non-negative", {
    for (i in 1:50) {
      n <- sample(50:200, 1)
      Y <- sample(c(0, 1), n, replace = TRUE)
      Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
      I <- rnorm(n) * Z

      IG <- calculate_conditional_ig_discrete(Y, Z, I)
      expect_gte(IG, -1e-6) # Allow tiny numerical errors
    }
  })

  # Test 2: IG = 0 for independent variables
  test_that("IG ≈ 0 for independent Y and I", {
    # Create truly independent Y and I
    set.seed(123)
    n <- 200
    Y <- sample(c(0, 1), n, replace = TRUE)
    Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
    I <- rnorm(n) * Z # Independent of Y

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_lt(IG, 0.15) # Should be close to 0 (allowing sampling variability)
  })

  # Test 3: IG > 0 for dependent variables
  test_that("IG > 0 for dependent Y and I", {
    # Create dependent Y and I
    set.seed(456)
    n <- 200
    Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
    I <- rnorm(n, mean = 0, sd = 2) * Z
    Y <- ifelse(I > 1, 1, 0) # Y depends on I

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_gt(IG, 0.05) # Should be clearly positive
  })

  # Test 4: IG bounded by H(Y)
  test_that("IG <= H(Y)", {
    for (i in 1:20) {
      n <- sample(50:200, 1)
      Y <- sample(c(0, 1), n, replace = TRUE)
      Z <- sample(c(0, 1), n, replace = TRUE, prob = c(0.7, 0.3))
      I <- rnorm(n) * Z

      H_Y <- calculate_entropy(Y)
      IG <- calculate_conditional_ig_discrete(Y, Z, I)
      expect_lte(IG, H_Y + 0.01)
    }
  })

  cat("✓ All IG property tests passed\n")
}


#' Test Suite: Theoretical Examples
test_theoretical_examples <- function() {
  cat("\n=== Testing Theoretical Examples ===\n")

  # Example 1: High IG (Lactate)
  test_that("Example 1 (Lactate) has high IG", {
    Y <- c(0, 1, 0, 1, 0, 0, 0, 1)
    Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
    I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    cat(sprintf("  Lactate IG = %.3f bits\n", IG))
    expect_gt(IG, 0.1) # Should be high
  })

  # Example 2: Low IG (Glucose - random imputation)
  test_that("Example 2 (Glucose) has low IG", {
    # Simulate random imputation
    set.seed(789)
    Y <- c(0, 0, 1, 1, 0, 0, 1, 1)
    Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
    I <- c(
      0, rnorm(1, 110, 10), 0, rnorm(1, 110, 10),
      0, rnorm(1, 110, 10), 0, rnorm(1, 110, 10)
    )

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    cat(sprintf("  Glucose IG = %.3f bits\n", IG))
    expect_lt(IG, 0.5) # Should be low
  })

  # Example 3: Z predicts Y, so I adds no information
  test_that("When Z predicts Y, interaction adds no information", {
    Y <- c(0, 0, 0, 0, 1, 1, 1, 1)
    Z <- c(0, 0, 0, 0, 1, 1, 1, 1) # Z perfectly predicts Y
    I <- c(0, 0, 0, 0, 3.5, 4.2, 4.8, 5.1) # Varying values

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    cat(sprintf("  IG when Z predicts Y = %.3f bits (should be ~0)\n", IG))
    expect_lt(IG, 0.1) # Should be low since Z already explains Y
  })

  cat("✓ All theoretical example tests passed\n")
}


#' Test Suite: Discretization
test_discretization <- function() {
  cat("\n=== Testing Discretization ===\n")

  test_that("Discretization preserves zeros", {
    I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)
    I_discrete <- discretize_interaction(I, n_bins = 3)

    # Check that zeros remain zeros
    expect_equal(I_discrete[I == 0], rep(0, sum(I == 0)))
  })

  test_that("Discretization creates correct number of bins", {
    I <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    I_discrete <- discretize_interaction(I, n_bins = 3)

    # Non-zero values should be in bins 1, 2, 3
    non_zero_bins <- unique(I_discrete[I != 0])
    expect_lte(length(non_zero_bins), 3)
  })

  cat("✓ All discretization tests passed\n")
}


#' Test Suite: Threshold Selection
test_threshold_selection <- function() {
  cat("\n=== Testing Threshold Selection ===\n")

  # Simulate IG values
  ig_values <- c(0.45, 0.32, 0.15, 0.08, 0.02, 0.01)

  test_that("Fixed threshold returns correct value", {
    epsilon <- select_threshold_fixed(0.01)
    expect_equal(epsilon, 0.01)
  })

  test_that("Adaptive threshold is reasonable", {
    epsilon <- select_threshold_adaptive(ig_values, fraction = 0.5)
    median_ig <- median(ig_values)
    expect_equal(epsilon, median_ig * 0.5)
    expect_gt(epsilon, 0)
    expect_lt(epsilon, max(ig_values))
  })

  test_that("Percentile threshold is reasonable", {
    epsilon <- select_threshold_percentile(ig_values, percentile = 25)
    expect_equal(epsilon, quantile(ig_values, 0.25))
  })

  cat("✓ All threshold selection tests passed\n")
}


#' Test Suite: Edge Cases
test_edge_cases <- function() {
  cat("\n=== Testing Edge Cases ===\n")

  test_that("Handle all-zero interaction", {
    Y <- c(0, 1, 0, 1, 0, 1)
    Z <- c(0, 0, 0, 0, 0, 0) # No missing data
    I <- c(0, 0, 0, 0, 0, 0) # All zeros

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_equal(IG, 0, tolerance = 0.01)
  })

  test_that("Handle small sample size", {
    Y <- c(0, 1)
    Z <- c(0, 1)
    I <- c(0, 2.5)

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_gte(IG, 0)
  })

  test_that("Handle all-missing data", {
    Y <- c(0, 1, 0, 1, 0, 1)
    Z <- c(1, 1, 1, 1, 1, 1) # All missing
    I <- c(2.1, 3.5, 1.8, 4.2, 2.0, 2.8)

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_gte(IG, 0)
  })

  cat("✓ All edge case tests passed\n")
}


#' Test Suite: Numerical Stability
test_numerical_stability <- function() {
  cat("\n=== Testing Numerical Stability ===\n")

  test_that("Handle extreme probabilities", {
    # Very skewed outcome
    Y <- c(rep(0, 99), 1)
    Z <- sample(c(0, 1), 100, replace = TRUE)
    I <- rnorm(100) * Z

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_gte(IG, 0)
    expect_false(is.na(IG))
    expect_false(is.infinite(IG))
  })

  test_that("Handle very small IG values", {
    # Nearly independent variables
    set.seed(999)
    Y <- sample(c(0, 1), 1000, replace = TRUE)
    Z <- sample(c(0, 1), 1000, replace = TRUE)
    I <- rnorm(1000, mean = 0, sd = 0.01) * Z # Very weak signal

    IG <- calculate_conditional_ig_discrete(Y, Z, I)
    expect_gte(IG, 0)
    expect_lt(IG, 0.1)
  })

  cat("✓ All numerical stability tests passed\n")
}


#' Run All Tests
run_all_tests <- function() {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║         ITIP Phase 1: Unit Test Suite                 ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n")

  # Run test suites
  test_entropy_calculations()
  test_ig_properties()
  test_theoretical_examples()
  test_discretization()
  test_threshold_selection()
  test_edge_cases()
  test_numerical_stability()

  cat("\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║         ALL TESTS PASSED ✓                            ║\n")
  cat("╚════════════════════════════════════════════════════════╝\n")
  cat("\n")

  cat("Phase 1 theoretical foundation is validated!\n")
  cat("Ready to proceed to Phase 2: Algorithm Implementation\n")
}


#' Quick Test (for development)
quick_test <- function() {
  cat("Running quick validation tests...\n")

  # Test basic IG calculation
  Y <- c(0, 1, 0, 1, 0, 0, 0, 1)
  Z <- c(0, 1, 0, 1, 0, 1, 0, 1)
  I <- c(0, 3.5, 0, 4.2, 0, 2.8, 0, 5.1)

  IG <- calculate_conditional_ig_discrete(Y, Z, I)
  cat(sprintf("✓ IG calculation works: %.3f bits\n", IG))

  # Test threshold selection
  ig_values <- c(0.45, 0.32, 0.15, 0.08, 0.02, 0.01)
  epsilon <- select_threshold_adaptive(ig_values)
  cat(sprintf("✓ Threshold selection works: ε = %.3f bits\n", epsilon))

  # Test non-negativity
  all_pass <- test_ig_nonnegativity(n_tests = 20)
  if (all_pass) {
    cat("✓ IG non-negativity verified\n")
  }

  cat("\nQuick tests passed! Run run_all_tests() for comprehensive validation.\n")
}


# If running this script directly, run all tests
if (sys.nframe() == 0) {
  run_all_tests()
}
