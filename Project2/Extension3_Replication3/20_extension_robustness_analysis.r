
# ==================================================
# Extension 3: Robustness analysis
# Step 1. Load and summarize independent network pools
# ==================================================

library(tidyverse)

# Prefer saved combined robustness results. When absent, combine
# the 100 original checkpoints and save a reusable copy.
combined_robustness_file <-
  "results/extension3/robustness_results_full.rds"

if (file.exists(combined_robustness_file)) {
  robustness_results <- readRDS(combined_robustness_file)
} else {
  robustness_files <- list.files(
    "results/extension3/robustness_simulation",
    pattern = "^chunk_[0-9]+\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  # Four new pools x five attachment levels x five checkpoints.
  stopifnot(length(robustness_files) == 4L * 5L * 5L)
  robustness_results <- purrr::map_dfr(robustness_files, readRDS)
  saveRDS(robustness_results, combined_robustness_file)
}

# Each of the four new pools has 500 Monte Carlo repetitions per
# attachment level and estimator. Check counts and unique rep IDs.
new_counts <- robustness_results %>%
  group_by(pool, attachment, estimator) %>%
  summarise(
    n_reps = n(),
    distinct_reps = n_distinct(rep),
    .groups = "drop"
  )

stopifnot(
  nrow(new_counts) == 4L * 5L * 3L,
  all(new_counts$n_reps == 500L),
  all(new_counts$distinct_reps == 500L)
)

# Load original main simulation results
extension_results <- readRDS(
  "results/extension3/extension_results_full.rds"
)

# Extract matching adoption setting from original pools
original_results <- extension_results %>%
  filter(
    alpha == -2,
    beta == 2,
    gamma == 0
  ) %>%
  mutate(pool = "pool_1")

# The original pool has 2,000 repetitions per estimator and
# attachment level; each additional pool has 500.
original_counts <- original_results %>%
  group_by(pool, attachment, estimator) %>%
  summarise(
    n_reps = n(),
    distinct_reps = n_distinct(rep),
    .groups = "drop"
  )
stopifnot(
  nrow(original_counts) == 5L * 3L,
  all(original_counts$n_reps == 2000L),
  all(original_counts$distinct_reps == 2000L)
)

# Combine original and new pools
robustness_all <- bind_rows(
  original_results,
  robustness_results
)

# Calculate estimator variances for each pool

# ==================================================
# Corrected robustness summary
# ==================================================

robustness_summary <- robustness_all %>%
  group_by(pool, attachment, estimator) %>%
  summarise(
    n_reps = n(),
    empirical_variance = var(estimate),
    .groups = "drop"
  ) %>%
  select(
    pool, attachment, estimator, empirical_variance
  ) %>%
  pivot_wider(
    id_cols = c(pool, attachment),
    names_from = estimator,
    values_from = empirical_variance
  ) %>%
  mutate(
    RE_Hajek = DM / Hajek,
    RE_HT = DM / HT
  )

print(
  robustness_summary %>%
    select(pool, attachment, RE_Hajek, RE_HT),
  n = Inf
)

# Save results used for the independent-pool robustness discussion.
output_file <- "results/extension3/robustness_summary.csv"
write_csv(robustness_summary, output_file)
cat("Saved robustness summary to:", output_file, "\n")




