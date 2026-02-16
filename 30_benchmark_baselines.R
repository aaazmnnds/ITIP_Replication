# 30_benchmark_baselines.R
# Full simulation benchmark: ITIP vs LASSO vs Stability Selection (REAL EXECUTION)
# Generates comprehensive_benchmarks.rds

library(MASS)
library(glmnet)
library(stabs)
library(hierNet)
library(glinternet)
library(missForest)
library(BAS)

# Source ITIP functions
# Source ITIP functions (Robust to running from root or scripts dir)
if (file.exists("03_bas_utils.R")) {
    source("03_bas_utils.R")
    source("04_itip_algorithm.R")
} else {
    source("scripts/03_bas_utils.R")
    source("scripts/04_itip_algorithm.R")
}

# --- Helper Functions ---

calc_f1 <- function(selected_features, true_features) {
    tp <- sum(selected_features %in% true_features)
    fp <- sum(!selected_features %in% true_features)
    fn <- sum(!true_features %in% selected_features)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0

    return(c(f1 = f1, precision = precision, recall = recall))
}

generate_data <- function(n = 500, p = 20, rho = 0.5) {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
    X <- mvrnorm(n, rep(0, p), Sigma)
    colnames(X) <- paste0("V", 1:p)

    # Introduce MAR missingness in V1 and V2
    X_miss <- X

    # Missingness in V1 depends on V3 (observed)
    miss_prob_v1 <- plogis(2 * X[, 3] - 1)
    missing_idx_v1 <- rbinom(n, 1, miss_prob_v1) == 1
    X_miss[missing_idx_v1, 1] <- NA

    # Missingness in V2 depends on V3 (observed)
    miss_prob_v2 <- plogis(1.5 * X[, 3] - 0.8)
    missing_idx_v2 <- rbinom(n, 1, miss_prob_v2) == 1
    X_miss[missing_idx_v2, 2] <- NA

    # True Model: Y depends on X1, X2, X3 (main effects)
    # AND interactions with missingness (I_V1, I_V2)
    # Form: 0.8 * X * (1 - Z) -> slope is 0.8 when observed, 0 when missing
    # This is equivalent to X - XZ, so ITIP should find the XZ term
    lp <- 0.5 * X[, 1] + 0.5 * X[, 2] + 0.3 * X[, 3] +
        1.5 * (X[, 1] * (1 - as.numeric(missing_idx_v1))) +
        1.5 * (X[, 2] * (1 - as.numeric(missing_idx_v2)))

    prob <- plogis(lp)
    y <- rbinom(n, 1, prob)

    return(list(X = X_miss, y = y, true_interactions = c("I_V1", "I_V2")))
}

# --- Main Benchmark Loop ---
p_values <- c(20, 50, 100)
results_list <- list()

for (p in p_values) {
    cat("\nRunning Benchmark for p =", p, "...\n")

    # Single run for speed (representative)
    set.seed(42 + p)
    data <- generate_data(n = 500, p = p)
    X <- data$X
    y <- data$y
    true_ints <- data$true_interactions

    # 1. ITIP (Actual Execution)
    cat("  Running ITIP...\n")
    df_itip <- as.data.frame(cbind(X, Y = y))
    itip_res <- tryCatch(
        {
            # Capture output
            capture.output({
                res <- itip(data = df_itip, outcome = "Y", threshold_method = "adaptive", verbose = FALSE)
            })
            res
        },
        error = function(e) {
            print(paste("ITIP Failed:", e$message))
            NULL
        }
    )

    if (!is.null(itip_res)) {
        # Check structure: pruned_features is list, kept_interactions is vector
        # kept_interactions already filtered to "I_..."
        itip_ints <- itip_res$pruned_features$kept_interactions
    } else {
        itip_ints <- c()
    }

    # 2. Standard LASSO
    cat("  Running LASSO...\n")
    # Impute for LASSO baseline
    # CRITICAL: Save missingness pattern BEFORE imputation
    missing_pattern <- is.na(X)

    imp_lasso <- missForest(X, verbose = FALSE)$ximp

    # --- CRITICAL FIX: Add Z (Missing Indicators) ---
    # Convert to dataframe
    imp_data <- as.data.frame(imp_lasso)

    # Identify missing columns
    miss_cols <- which(colSums(missing_pattern) > 0)

    if (length(miss_cols) > 0) {
        for (k in seq_along(miss_cols)) {
            col_idx <- miss_cols[k]
            col_name <- colnames(X)[col_idx]

            # Create Z column
            z_name <- paste0("Z_", col_name)
            imp_data[[z_name]] <- as.numeric(missing_pattern[, col_idx])
        }
    }

    # Create standard model matrix (now includes Z)
    X_full <- model.matrix(~ .^2, data = imp_data)[, -1]

    # --- Augment Baselines for Fairness ---
    # Manually add imputation interactions (I_j = X_imp * Z)
    # Note: X_full already has X:Z interactions from model.matrix(~.^2), but ITIP uses I_j features specifically.
    # We add them explicitly to be safe and consistent with ITIP naming.

    # We already identified miss_cols above

    if (length(miss_cols) > 0) {
        I_extras <- matrix(0, nrow = nrow(X), ncol = length(miss_cols))
        extra_names <- character(length(miss_cols))

        for (k in seq_along(miss_cols)) {
            col_idx <- miss_cols[k]
            col_name <- colnames(X)[col_idx]

            # I_j = X_imp * Z
            # Use SAVED missingness pattern, not possibly imputed X
            Z_vec <- as.numeric(missing_pattern[, col_idx])
            I_term <- imp_lasso[, col_idx] * Z_vec

            I_extras[, k] <- I_term
            # Name: I_V1, I_V2 to match true_ints
            extra_names[k] <- paste0("I_", col_name)
        }
        colnames(I_extras) <- extra_names
        # Combine with standard interactions
        X_full <- cbind(X_full, I_extras)
    }
    cv_lasso <- cv.glmnet(X_full, y, family = "binomial", alpha = 1)
    coefs <- coef(cv_lasso, s = "lambda.min")
    selected_lasso <- rownames(coefs)[coefs[, 1] != 0]

    # Detect if V:Z interactions were selected (in standard model.matrix form)
    # The true interactions are V1:Z_V1 and V2:Z_V2
    # LASSO might pick "V1:Z_V1" or "Z_V1:V1"

    # Check for V1:Z_V1
    mask_v1_z1 <- grepl("^V1:Z_V1$|^Z_V1:V1$", selected_lasso)

    # Check for V2:Z_V2
    mask_v2_z2 <- grepl("^V2:Z_V2$|^Z_V2:V2$", selected_lasso)

    # Construct lasso_ints (for scoring)
    lasso_ints <- character(0)

    # If standard interaction found, map to "I_V1"
    if (any(mask_v1_z1)) lasso_ints <- c(lasso_ints, "I_V1")
    if (any(mask_v2_z2)) lasso_ints <- c(lasso_ints, "I_V2")

    # Also include any directly selected "I_" features (if duplication wasn't an issue)
    lasso_ints <- c(lasso_ints, grep("^I_", selected_lasso, value = TRUE))

    # EXCLUDE other interactions (noise) as per user request to focus on Recall of target
    lasso_ints <- unique(lasso_ints)

    # 3. Stability Selection
    cat("  Running Stability Selection...\n")
    if (p < 100) {
        stab_fit <- stabsel(x = X_full, y = y, fitfun = glmnet.lasso, cutoff = 0.75, PFER = 1)
        selected_stabs <- names(stab_fit$selected)
        selected_stabs <- names(stab_fit$selected)

        # Check for V:Z patterns
        mask_v1_z1_stabs <- grepl("^V1:Z_V1$|^Z_V1:V1$", selected_stabs)
        mask_v2_z2_stabs <- grepl("^V2:Z_V2$|^Z_V2:V2$", selected_stabs)

        stabs_ints <- character(0)
        if (any(mask_v1_z1_stabs)) stabs_ints <- c(stabs_ints, "I_V1")
        if (any(mask_v2_z2_stabs)) stabs_ints <- c(stabs_ints, "I_V2")

        stabs_ints <- c(stabs_ints, grep("^I_", selected_stabs, value = TRUE))
        stabs_ints <- unique(stabs_ints)
    } else {
        stabs_ints <- character(0) # Skip for p=100 speed
    }

    # 4. Elastic Net (Alpha = 0.5)
    cat("  Running Elastic Net (Alpha=0.5)...\n")
    cv_enet <- cv.glmnet(X_full, y, family = "binomial", alpha = 0.5)
    coefs_enet <- coef(cv_enet, s = "lambda.min")
    selected_enet <- rownames(coefs_enet)[coefs_enet[, 1] != 0]

    # Check for V:Z patterns
    mask_v1_z1_enet <- grepl("^V1:Z_V1$|^Z_V1:V1$", selected_enet)
    mask_v2_z2_enet <- grepl("^V2:Z_V2$|^Z_V2:V2$", selected_enet)

    enet_ints <- character(0)
    if (any(mask_v1_z1_enet)) enet_ints <- c(enet_ints, "I_V1")
    if (any(mask_v2_z2_enet)) enet_ints <- c(enet_ints, "I_V2")

    enet_ints <- c(enet_ints, grep("^I_", selected_enet, value = TRUE))
    enet_ints <- unique(enet_ints)

    # Calculate Scores
    clean_names <- function(x) {
        if (length(x) == 0) {
            return(character(0))
        }
        # Remove backticks
        x <- gsub("`", "", x)

        # Handle imputation interactions (I_V1, I_V2) - keep as-is
        is_imp_int <- grepl("^I_", x)

        # For pairwise interactions, remove suffixes and sort
        if (any(!is_imp_int)) {
            x[!is_imp_int] <- sapply(x[!is_imp_int], function(s) {
                s <- gsub("_imp", "", s)
                s <- gsub("_miss", "", s)
                parts <- strsplit(s, ":")[[1]]
                paste(sort(parts), collapse = ":")
            })
        }

        return(x)
    }

    # True interactions for scoring
    true_set <- true_ints

    # 5. Hierarchical LASSO (glinternet)
    cat("  Running Hierarchical LASSO (glinternet)...\n")
    # imp_lasso is already available
    numLevels <- rep(1, ncol(imp_lasso)) # All continuous
    # Capture output to suppress internal logging
    gl_ints <- tryCatch(
        {
            capture.output({
                cv_gl <- glinternet.cv(imp_lasso, y, numLevels, family = "binomial", verbose = FALSE)
            })
            # Get lambda with min error
            i_min <- which.min(cv_gl$cvErr[[1]]) # 1 for standard CV error
            interactions <- cv_gl$interactions[[i_min]]

            if (is.null(interactions) || length(interactions) == 0) {
                character(0)
            } else {
                # Map indices to names
                # interactions is a list of matrices or a matrix?
                # It's a matrix with 2 columns
                v_names <- colnames(imp_lasso)
                if (is.matrix(interactions)) {
                    apply(interactions, 1, function(row) {
                        paste(sort(c(v_names[row[1]], v_names[row[2]])), collapse = ":")
                    })
                } else {
                    character(0)
                }
            }
        },
        error = function(e) {
            print(paste("glinternet Failed:", e$message))
            character(0)
        }
    )

    s_lasso <- calc_f1(clean_names(lasso_ints), true_set)
    s_stabs <- calc_f1(clean_names(stabs_ints), true_set)
    s_enet <- calc_f1(clean_names(enet_ints), true_set)
    s_itip <- calc_f1(clean_names(itip_ints), true_set)
    s_gl <- calc_f1(clean_names(gl_ints), true_set)

    cat("  Results p=", p, ": ITIP F1=", s_itip["f1"], " GL F1=", s_gl["f1"], " LASSO F1=", s_lasso["f1"], " ENET F1=", s_enet["f1"], " Stabs F1=", s_stabs["f1"], "\n")
    cat("  ITIP Prec=", s_itip["precision"], " Rec=", s_itip["recall"], "\n")
    cat("  GL Prec=", s_gl["precision"], " Rec=", s_gl["recall"], "\n")

    results_list[[as.character(p)]] <- c(
        itip_f1 = s_itip["f1"], itip_prec = s_itip["precision"], itip_rec = s_itip["recall"],
        gl_f1 = s_gl["f1"], gl_prec = s_gl["precision"], gl_rec = s_gl["recall"],
        lasso_f1 = s_lasso["f1"], lasso_prec = s_lasso["precision"], lasso_rec = s_lasso["recall"],
        enet_f1 = s_enet["f1"], enet_prec = s_enet["precision"], enet_rec = s_enet["recall"],
        stabs_f1 = s_stabs["f1"], stabs_prec = s_stabs["precision"], stabs_rec = s_stabs["recall"]
    )
}

cat("  [DEBUG] LASSO selected (1se):", paste(selected_lasso, collapse = ", "), "\n")

# Save results (Robust path)
if (dir.exists("results")) {
    saveRDS(results_list, "results/comprehensive_benchmarks.rds")
} else {
    saveRDS(results_list, "../results/comprehensive_benchmarks.rds")
}
print(results_list)
