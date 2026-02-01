# ITIP Simulation Diagnostic

# Goal: Test if conditioning on X_imp fixes the "all significant" problem in v3

library(BAS)
library(missForest)
library(MASS)

# 1. Setup (Same as v3)
set.seed(42)
n <- 2000
p <- 5
beta_int <- 4.0
missing_rate <- 0.3
rho <- 0.6

# Generate Data
mu <- rep(0, p)
Sigma <- matrix(rho, nrow = p, ncol = p)
diag(Sigma) <- 1.0
X_full <- mvrnorm(n, mu, Sigma)
colnames(X_full) <- paste0("X", 1:p)

Z <- matrix(0, nrow = n, ncol = p)
for (j in 1:p) Z[, j] <- rbinom(n, 1, missing_rate)

X_incomp <- X_full
X_incomp[Z == 1] <- NA
imp_res <- missForest(X_incomp, verbose = FALSE, ntree = 50) # Faster
X_imp <- imp_res$ximp

I <- matrix(0, nrow = n, ncol = p)
for (j in 1:p) I[, j] <- X_imp[, j] * Z[, j]

# Outcome (True I_3)
logits_y <- 0.5 * X_full[, 1] - 0.5 * X_full[, 2] + beta_int * (X_full[, 3] * Z[, 3])
probs_y <- 1 / (1 + exp(-logits_y))
Y <- rbinom(n, 1, probs_y)

# 2. Test Corrected Model: Y ~ X_imp + Z + I vs Y ~ X_imp + Z
cat("\nDiagnostic Results (Conditioning on X_imp):\n")
for (j in 1:p) {
    # Full Model: X_imp + Z + I
    design <- data.frame(Y = Y, X_j = X_imp[, j], Z_j = Z[, j], I_j = I[, j])

    # We want to see if I_j is needed given X_j and Z_j
    bas_fit <- bas.glm(
        Y ~ X_j + Z_j + I_j,
        data = design, family = binomial(),
        method = "MCMC", n.models = 2000, MCMC.iterations = 5000,
        include.always = ~ X_j + Z_j
    )

    # Check PIP of I_j (4th term)
    pip <- bas_fit$probne0[4] # Intercept, X, Z, I
    ig <- -log(1 - pip + 1e-6)

    cat(sprintf("I_%d: PIP=%.3f, IG=%.3f\n", j, pip, ig))
}
