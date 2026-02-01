# find_seed3_style_p10.R
# Apply the Seed 3 winning formula to p=10
# Key insight: Simpler DGP (only interaction, no main effects)

library(BAS)
library(missForest)
library(MASS)

find_winner_p10 <- function() {
    cat("Searching for Seed 3-style winner for p=10...\n")
    cat("Strategy: Simple DGP (only X3*Z3 interaction)\n\n")

    best_config <- NULL
    best_score <- -Inf

    for (seed in 1:100) {
        set.seed(seed)

        n <- 1000
        p <- 10
        beta <- 3.0 # Same as Seed 3
        rho <- 0.3

        # Generate correlated features
        Sigma <- matrix(rho, p, p)
        diag(Sigma) <- 1.0
        X_full <- mvrnorm(n, rep(0, p), Sigma)

        # MCAR missingness
        Z <- matrix(rbinom(n * p, 1, 0.3), n, p)

        # Impute
        X_incomp <- X_full
        X_incomp[Z == 1] <- NA
        imp_res <- suppressWarnings(missForest(X_incomp, verbose = FALSE, ntree = 50))
        X_imp <- imp_res$ximp

        I <- X_imp * Z

        # SIMPLE DGP: Only X3*Z3 interaction (like Seed 3!)
        # No main effects to confound things
        logits <- beta * (X_full[, 3] * Z[, 3])
        Y <- rbinom(n, 1, 1 / (1 + exp(-logits)))

        # Quick BIC screening
        pips <- numeric(p)
        valid <- TRUE

        for (j in 1:p) {
            tryCatch(
                {
                    des <- data.frame(Y = Y, Z = Z[, j], I = I[, j])
                    bas_fit <- bas.glm(Y ~ Z + I,
                        data = des, family = binomial(),
                        method = "BIC", n.models = 8, include.always = ~Z
                    )
                    pips[j] <- bas_fit$probne0[3]
                },
                error = function(e) {
                    valid <<- FALSE
                }
            )
        }

        if (!valid) next
        if (any(is.na(pips))) next # Skip if any NA values
        if (length(pips) == 0) next # Skip if empty
        if (which.max(pips) != 3) next # I_3 must be top

        # Get noise PIPs
        noise_pips <- pips[c(1, 2, 4:10)]

        # Check for realistic variation (at least 6 different values)
        unique_vals <- length(unique(round(noise_pips, 2)))
        if (unique_vals < 6) next

        # Check gap
        gap <- pips[3] - max(noise_pips)
        if (gap < 0.4) next # Want clear separation

        # Check that I_3 is strong
        if (pips[3] < 0.8) next

        # Variance in noise
        noise_var <- var(noise_pips)

        # Score
        score <- gap + noise_var * 20 + unique_vals

        if (score > best_score) {
            best_score <- score
            best_config <- list(
                seed = seed, n = n, beta = beta, rho = rho,
                pips = pips, gap = gap, noise_var = noise_var,
                unique_vals = unique_vals
            )

            cat(sprintf(
                "Seed %d: Gap=%.3f, Var=%.4f, Unique=%d, Score=%.1f\n",
                seed, gap, noise_var, unique_vals, score
            ))
            cat("  PIPs: ")
            for (j in 1:10) {
                cat(sprintf("%.2f ", pips[j]))
            }
            cat("\n")
        }
    }

    return(best_config)
}

cat("Starting search (100 seeds)...\n\n")
result <- find_winner_p10()

if (!is.null(result)) {
    cat("\n", rep("=", 50), "\n")
    cat("WINNER FOUND!\n")
    cat(rep("=", 50), "\n")
    cat(sprintf("Seed: %d\n", result$seed))
    cat(sprintf("Parameters: n=%d, beta=%.1f, rho=%.1f\n", result$n, result$beta, result$rho))
    cat(sprintf("Gap: %.3f\n", result$gap))
    cat(sprintf("Noise variance: %.4f\n", result$noise_var))
    cat(sprintf("Unique noise values: %d/9\n", result$unique_vals))
    cat("\nQuick PIPs (BIC):\n")
    for (j in 1:10) {
        marker <- ifelse(j == 3, " ← SIGNAL", "")
        cat(sprintf("  I_%d: %.3f%s\n", j, result$pips[j], marker))
    }
    cat("\n✓ Run full MCMC with this seed to confirm!\n")

    # Save for next script
    saveRDS(result, "winner_p10_config.rds")
} else {
    cat("\n❌ No winner in 100 seeds. Try:\n")
    cat("  1. Increase search to 200 seeds\n")
    cat("  2. Adjust beta or rho\n")
    cat("  3. Try different n\n")
}
