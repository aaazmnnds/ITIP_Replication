# bas_utils.R
# Utility functions for Bayesian Adaptive Sampling (BAS) integration
# Author: Azman Nads

library(BAS)

#' Run BAS for Variable Selection
#'
#' Estimates Posterior Inclusion Probabilities (PIP) for features using
#' Bayesian Model Averaging.
#'
#' @param formula Model formula (e.g., Y ~ .)
#' @param data Data frame containing outcome and features
#' @param method Sampling method ("MCMC", "BAS", "deterministic")
#' @param n_models Number of models to sample (for MCMC/BAS)
#' @param prior Prior on coefficients ("BIC", "AIC", "g-prior", "hyper-g")
#' @param model_prior Prior on model space ("uniform", "beta-binomial")
#' @param verbose Print progress
#' @return List containing:
#'   - bas_object: Full BAS result object
#'   - pip: Named vector of Posterior Inclusion Probabilities
#'   - best_model: Features in the highest probability model
run_bas_selection <- function(formula, data,
                              method = "MCMC",
                              n_models = 10000,
                              prior = "hyper-g",
                              model_prior = "beta-binomial",
                              include.always = NULL,
                              verbose = FALSE) {
    if (verbose) cat("  Running BAS selection...\n")

    # Ensure data is clean (no missing values allowed in BAS input)
    if (any(is.na(data))) {
        stop("Data contains missing values. BAS requires complete/imputed data.")
    }

    # Run BAS
    # family = binomial() for binary outcomes (mortality/sepsis)
    # Construct betaprior object
    if (prior == "AIC") {
        prior_obj <- AIC()
    } else if (prior == "BIC") {
        prior_obj <- BIC()
    } else if (prior == "hyper-g") {
        prior_obj <- hyper.g.n() # Default hyper-g/n
    } else if (prior == "g-prior") {
        prior_obj <- g.prior()
    } else if (prior == "zellner-siow") {
        prior_obj <- JZS()
    } else {
        # Fallback or allow passing object directly if implemented later
        warning("Unknown prior string. Using hyper-g(n) as default.")
        prior_obj <- hyper.g.n()
    }

    bas_fit <- bas.glm(
        formula = formula,
        data = data,
        family = binomial(),
        method = method,
        n.models = n_models,
        betaprior = prior_obj,
        modelprior = beta.binomial(1, 1), # Default beta-binomial(1,1) is uniform on model size
        include.always = include.always,
        initprobs = "eplogp" # Initialize using p-values
    )

    # Extract PIPs
    # probne0 gives the marginal posterior inclusion probability
    pip <- bas_fit$probne0
    names(pip) <- bas_fit$namesx

    # Remove intercept from PIP list if present (usually index 1 but safe to check name)
    if ("Intercept" %in% names(pip)) {
        pip <- pip[names(pip) != "Intercept"]
    }

    # Identify best model (Highest Posterior MArginall Probability - HPM or MPM)
    # Here we usually look at MPM (Median Probability Model) -> PIP > 0.5
    mpm_features <- names(pip)[pip > 0.5]

    if (verbose) {
        cat(sprintf("  BAS identified %d features with PIP > 0.5\n", length(mpm_features)))
    }

    return(list(
        bas_object = bas_fit,
        pip = pip,
        mpm_features = mpm_features
    ))
}
