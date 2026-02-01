# check_benchmarks.R
# Read RDS files and print summary for p=20, 50, 100

results <- tryCatch(
    {
        readRDS("../results/comprehensive_benchmarks.rds")
    },
    error = function(e) {
        NULL
    }
)

results_large <- tryCatch(
    {
        readRDS("../results/large_scale_benchmarks.rds")
    },
    error = function(e) {
        NULL
    }
)

print("--- Standard Results (p=20, 50, 100) ---")
if (!is.null(results)) print(results) else print("Standard RDS not found/readable.")

print("--- Large Scale Results ---")
if (!is.null(results_large)) print(results_large) else print("Large Scale RDS not found/readable.")
