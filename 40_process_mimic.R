# process_mimic.R
# Complete MIMIC-III Septic Shock Pipeline: Merge -> Clean -> Finalize
# Author: Azman Nads

library(dplyr)

# === Configuration ===
# Input Paths
tabular_path <- "../data/mimic_septic_shock_tabular.csv"
admissions_path <- "../data/ADMISSIONS.csv.gz"

# Output Path
output_path <- "../data/mimic_septic_shock_final.csv"

cat("=== MIMIC-III Data Processing Pipeline ===\n")

# 1. Validation
if (!file.exists(tabular_path)) stop(sprintf("Tabular file not found: %s", tabular_path))
if (!file.exists(admissions_path)) stop(sprintf("Admissions file not found: %s", admissions_path))

# 2. Load Raw Tabular Data
cat("Step 1: Loading clinical features...\n")
data <- read.csv(tabular_path, stringsAsFactors = FALSE)
cat(sprintf("  Loaded %d rows, %d cols\n", nrow(data), ncol(data)))

# 3. Load & Merge Outcome (Mortality)
cat("Step 2: Merging mortality labels from ADMISSIONS...\n")
df_admissions <- read.csv(gzfile(admissions_path), stringsAsFactors = FALSE)
df_admissions <- df_admissions[, c("HADM_ID", "HOSPITAL_EXPIRE_FLAG")]

# Merge
data_merged <- merge(data, df_admissions, by = "HADM_ID", all.x = TRUE)

# Rename outcome
names(data_merged)[names(data_merged) == "HOSPITAL_EXPIRE_FLAG"] <- "mortality"

# Drop rows with missing outcome
n_missing <- sum(is.na(data_merged$mortality))
if (n_missing > 0) {
    warning(sprintf("  Dropping %d rows with missing mortality labels", n_missing))
    data_merged <- data_merged[!is.na(data_merged$mortality), ]
}
cat(sprintf(
    "  Result: %d rows. Mortality Rate: %.2f%%\n",
    nrow(data_merged), 100 * mean(data_merged$mortality)
))

# 4. Cleaning & Feature Selection
cat("Step 3: Cleaning headers and features...\n")
# Drop IDs and target leaks
drop_cols <- c("SPECIMEN_TYPE", "DEATHTIME", "subject_id", "stay_id")
cols_to_drop <- intersect(names(data_merged), drop_cols)

if (length(cols_to_drop) > 0) {
    data_final <- data_merged[, !(names(data_merged) %in% cols_to_drop)]
    cat(sprintf("  Dropped operational columns: %s\n", paste(cols_to_drop, collapse = ", ")))
} else {
    data_final <- data_merged
}

# 5. Save Final Dataset
cat("Step 4: Saving final dataset...\n")
write.csv(data_final, output_path, row.names = FALSE)
cat(sprintf("  Saved to %s\n", output_path))
cat("=== Processing Complete ===\n")
