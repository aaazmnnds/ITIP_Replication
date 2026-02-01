library(BAS)
library(missForest)
library(ggplot2)

# MNAR Robustness Simulation
# Reviewer concerns: Can ITIP find interactions when Z itself is a strong predictor?

run_mnar_robustness <- function(n = 500, alpha_values = seq(0, 3, by = 0.5)) {
    results <- data.frame()

    for (alpha in alpha_values) {
        cat("Testing alpha =", alpha, "\n")

        # 1. Generate Data
        set.seed(42)
        x <- rnorm(n)
        # MNAR: Z depends on X
        prob_z <- plogis(x - 0.5)
        z <- rbinom(n, 1, prob_z)

        # Outcome Y: logit(Y) = 0.5*X + alpha*Z + 2.0*(X*Z)
        # alpha is the 'masking' signal of the indicator
        logits <- 0.5 * x + alpha * z + 2.0 * (x * z)
        y <- rbinom(n, 1, plogis(logits))

        # 2. ITIP Application
        x_obs <- x
        x_obs[z == 1] <- NA
        df <- data.frame(x = x_obs, y = y, dummy = rnorm(n))

        # Impute
        cat("  Missing values count:", sum(is.na(df$x)), "\n")
        cat("  Complete rows count:", sum(!is.na(df$x)), "\n")
        if (any(is.na(df)) && sum(!is.na(df$x)) > 5) {
            imp <- missForest(df, ntree = 10)$ximp
        } else {
            cat("  Bypassing missForest (no missing values or too few complete rows)\n")
            # If all missing, use mean imputation as simple fallback
            if (any(is.na(df))) {
                df$x[is.na(df$x)] <- mean(x, na.rm = TRUE)
            }
            imp <- df
        }

        # Construct Interaction
        int_df <- data.frame(y = y, x_imp = imp$x, z = z, i = imp$x * z)

        # BAS Estimation
        model <- bas.glm(y ~ z + i,
            data = int_df, family = binomial(), method = "BAS",
            update = 200, include.always = ~z
        )
        # We always include z (as per ITIP logic) and check pip for i
        pip_i <- model$probne0[3] # Index 3 is 'i'
        ig_proxy <- -log(1 - min(pip_i, 0.9999))

        results <- rbind(results, data.frame(alpha = alpha, pip_i = pip_i, ig_proxy = ig_proxy))
    }

    return(results)
}

# Run and Plot
mnar_results <- run_mnar_robustness()
saveRDS(mnar_results, "../results/mnar_robustness_results.rds")

# Plotting
p <- ggplot(mnar_results, aes(x = alpha, y = ig_proxy)) +
    geom_line() +
    geom_point() +
    labs(
        title = "ITIP Robustness to MNAR 'Masking'",
        subtitle = "Information Gain of Interaction (i) as Indicator Signal (alpha) Increases",
        x = "Strength of Missing Indicator Signal (alpha)",
        y = "ITIP IG Proxy (-log(1-gamma))"
    ) +
    theme_minimal()

ggsave("../results/mnar_robustness_plot.png", p)
print(mnar_results)
