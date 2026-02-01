# generate_simulation_data.R
# Generate synthetic datasets for ITIP validation (Simulation Study)
# Author: Azman Nads

#' Generate Synthetic Data for ITIP
#'
#' Creates a dataset with known ground-truth interactions and missingness mechanisms.
#'
#' @param n Sample size
#' @param p Number of features
#' @param n_interactions Number of true contributing interactions
#' @param missingness_rate Proportion of missing summaries (e.g. 0.3)
#' @param mechanism Missingness mechanism ("MCAR", "MNAR")
#' @param seed Random seed
#' @return List containing:
#'   - data: Data frame with missing values (features + Y)
#'   - data_complete: Complete data frame (before missingness)
#'   - truth: List of true features (important for validation)
generate_simulation_data <- function(n = 500, p = 10, n_interactions = 2,
                                     missingness_rate = 0.3, mechanism = "MNAR",
                                     seed = 42) {
    set.seed(seed)

    # 1. Generate Features (Correlated Gaussian)
    # Simple: Independent for now, or use MASS::mvrnorm for correlation
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- paste0("X", 1:p)

    # 2. Define True Model
    # Let's say Y depends on:
    # - X1 (Main effect)
    # - X2 (Main effect)
    # - Interaction: X3 * Z_X3 (where Z_X3 is the missing indicator for X3)
    #   This simulates a case where "missingness" itself modifies the effect of the value (imputed).
    #   Wait, in ITIP theory: Interaction is between Imputed Value and Missing Indicator.
    #   Let's simulate a simpler "biological" interaction first, then "missingness" interaction.

    # For ITIP validation, we specifically want to test if it finds:
    # Interaction: X_imp * Z

    # Let's force X3 to have missingness.
    # Interaction Effect: The outcome Y is higher if X3 is HIGH AND X3 is MISSING (observed via imputation? No, Z=1).
    # Actually, "Imputed Value * Indicator" is tricky to simulate in ground truth because Z doesn't exist yet.
    # But the generative process is:
    # Y depends on X3 (latent) and Z3 (missing status).
    # Y ~ X1 + X2 + beta * (X3 * Z3) -> This is the interaction we want ITIP to find.

    # Pre-define missingness pattern (Z)
    Z_matrix <- matrix(0, nrow = n, ncol = p)

    if (mechanism == "MCAR") {
        # Randomly remove values
        for (j in 1:p) {
            mis_indices <- sample(1:n, size = floor(missingness_rate * n))
            Z_matrix[mis_indices, j] <- 1
        }
    } else if (mechanism == "MNAR") {
        # Missing dependent on value (e.g., higher values missing)
        for (j in 1:p) {
            # Probabilistic missingness based on value
            prob_missing <- 1 / (1 + exp(-(X[, j] - 0.5))) # Logistic function
            Z_matrix[, j] <- rbinom(n, 1, prob_missing)
        }
    }

    # Generate Outcome Y
    # Model: Logit(P(Y=1)) = X1 + X2 + 2 * (X3 * Z3)
    # Meaning: X3 affects Y ONLY (or differently) when it is missing?
    # Or "If missing, risk is higher proportional to the underlying value"?
    # Let's simplify:
    # Main effects: X1, X2
    # Interaction: X3 * Z3 (The value of X3 matters, but only when it was missing?
    #   No, that implies we know X3 when it's missing. Standard imputation fills it.
    #   Let's assume "Z3" is a risk factor, and it interacts with the imputed value.)

    # True logits (using complete data X)
    # We use X3 * Z3. Since X3 is latent when Z3=1, this represents the "true" value * missingness interaction.
    logits <- 0.5 * X[, 1] - 0.5 * X[, 2] +
        2.0 * (X[, 3] * Z_matrix[, 3]) # Strong interaction target

    probs <- 1 / (1 + exp(-logits))
    Y <- rbinom(n, 1, probs)

    # 3. Create Observed Data (Mask values)
    X_obs <- X
    X_obs[Z_matrix == 1] <- NA

    # Return
    data <- as.data.frame(X_obs)
    data$Y <- Y

    data_complete <- as.data.frame(X)
    data_complete$Y <- Y

    truth <- list(
        main_effects = c("X1", "X2"),
        interactions = c("I_X3"), # I_X3 corresponds to X3 * Z3
        Z_cols = paste0("Z_X", 1:p)[colSums(Z_matrix) > 0]
    )

    return(list(
        data = data,
        data_complete = data_complete,
        truth = truth,
        Z_matrix = Z_matrix
    ))
}

# Example generation
if (FALSE) {
    sim <- generate_simulation_data(n = 200, missingness_rate = 0.3)
    print(head(sim$data))
}
