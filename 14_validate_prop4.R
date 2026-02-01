# Validation script for Proposition 4 (IG Proxy)
library(BAS)
library(ggplot2)

# Manual Mutual Information for discrete X, Y
mutinformation <- function(x, y) {
    p_xy <- table(x, y) / length(x)
    p_x <- rowSums(p_xy)
    p_y <- colSums(p_xy)

    mi <- 0
    for (i in seq_len(nrow(p_xy))) {
        for (j in seq_len(ncol(p_xy))) {
            if (p_xy[i, j] > 0) {
                mi <- mi + p_xy[i, j] * log(p_xy[i, j] / (p_x[i] * p_y[j]))
            }
        }
    }
    return(mi)
}

discretize <- function(x, bins = 10) {
    return(cut(x, breaks = bins, labels = FALSE))
}

set.seed(42)

# Parameters
n <- 1000
p <- 100
snrs <- seq(0, 1, length.out = p) # Smaller signal strengths to avoid saturation

# Data Generation
generate_data <- function(n, snr) {
    x <- rnorm(n)
    # Outcome Y depends on X with varying strength
    prob <- plogis(snr * x)
    y <- rbinom(n, 1, prob)
    return(data.frame(x = x, y = y))
}

results <- data.frame()

for (i in 1:p) {
    snr <- snrs[i]
    df <- generate_data(n, snr)

    # 1. Compute "True" Mutual Information (discretized)
    # discretize Y is already 0/1, discretize X
    x_disc <- discretize(df$x)
    true_mi <- mutinformation(x_disc, df$y)

    # 2. Compute BAS Posterior Probability
    # For BAS, we fit a model Y ~ X
    model <- bas.glm(y ~ x, data = df, family = binomial(), method = "BAS", update = 100)
    gamma_j <- model$probne0[2] # probability for 'x'

    # 3. Compute Proxy
    # Avoid log(0)
    proxy_val <- -log(1 - min(gamma_j, 0.9999))

    results <- rbind(results, data.frame(snr = snr, true_mi = true_mi, gamma_j = gamma_j, proxy = proxy_val))
}

# Plot
p_plot <- ggplot(results, aes(x = true_mi, y = proxy)) +
    geom_point(aes(color = snr), size = 3) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed") +
    labs(
        title = "Validation of Proposition 4: Proxy vs True MI",
        x = "Mutual Information (discretized)",
        y = "-log(1 - gamma_j) Proxy",
        color = "SNR"
    ) +
    theme_minimal()

print(p_plot)

# Correlation
cor_val <- cor(results$true_mi, results$proxy, method = "spearman")
cat("Spearman Correlation between True MI and Proxy:", cor_val, "\n")

# Save findings
write.csv(results, "../results/prop4_validation_results.csv", row.names = FALSE)
ggsave("../results/prop4_validation_plot.png", p_plot)
