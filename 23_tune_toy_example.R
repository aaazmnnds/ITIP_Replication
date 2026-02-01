# tune_toy_example.R
# Hypothesis: Increasing n make p=10 possible, provided we lower beta to prevent saturation.
# Search Space:
#   n: {2000, 3000}  <- Much larger sample size
#   beta: {1.0, 1.5, 2.0} <- Weaker signal (relying on n to find it)
#   rho: {0.3} <- Fixed moderate correlation

library(BAS)
library(missForest)
library(MASS)

find_large_n_winner <- function() {
    cat("Searching for p=10 configuration with LARGE N...\n")

    # Grid search
    n_list <- c(2000, 3000)
    beta_list <- c(1.0, 1.5, 2.0)

    best_config <- NULL
    best_score <- -Inf

    for (n in n_list) {
        for (beta in beta_list) {
            # Try 20 seeds per configuration
            for (seed in 1:20) {
                set.seed(seed)

                p <- 10
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

                # SIMPLE DGP: Only X3*Z3 interaction
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
                if (any(is.na(pips))) next
                if (length(pips) == 0) next

                # CRITERIA:
                # 1. I_3 is top
                if (which.max(pips) != 3) next

                # 2. Strong Signal
                if (pips[3] < 0.90) next

                # 3. Good Separation
                noise_pips <- pips[c(1, 2, 4:10)]
                max_noise <- max(noise_pips)
                gap <- pips[3] - max_noise

                if (gap < 0.5) next

                # 4. Realistic Noise (not all identical)
                unique_vals <- length(unique(round(noise_pips, 2)))
                if (unique_vals < 5) next

                score <- gap

                if (score > best_score) {
                    best_score <- score
                    best_config <- list(seed = seed, n = n, beta = beta, pips = pips, gap = gap)

                    cat(sprintf(
                        "Found Candidate: N=%d Beta=%.1f Seed=%d -> Gap=%.3f (I3=%.3f, MaxN=%.3f)\n",
                        n, beta, seed, gap, pips[3], max_noise
                    ))

                    # Stop if perfect
                    if (gap > 0.8 && unique_vals >= 6) {
                        return(best_config)
                    }
                }
            }
        }
    }
    return(best_config)
}

res <- find_large_n_winner()

if (!is.null(res)) {
    cat("\n=== WINNER PARAMETERS ===\n")
    cat(sprintf("Seed: %d\n", res$seed))
    cat(sprintf("N: %d\n", res$n))
    cat(sprintf("Beta: %.1f\n", res$beta))
    cat("\nPIPs:\n")
    for (j in 1:10) {
        marker <- ifelse(j == 3, " *", "")
        cat(sprintf("I_%d: %.3f%s\n", j, res$pips[j], marker))
    }
    saveRDS(res, "winner_p10_large_n.rds")
} else {
    cat("No configuration found.\n")
}
