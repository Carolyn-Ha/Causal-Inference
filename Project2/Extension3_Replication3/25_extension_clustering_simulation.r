
# ============================================================
# 25_extension_clustering_simulation.r
#
# Degree-preserving topology experiment
#
# Compare 150 matched low- and high-clustering networks.
#
# Analysis A:
#   Estimate treatment effects within each matched pair.
#
# Analysis B:
#   Compare DM, HT, and Hajek estimator variances
#   across repeated 150-village experiments.
#
# Uses our EXISTING adoption model and seeding functions.
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(tidyverse)
library(igraph)

# Run this script from the project root (where 01_...r through 25_...r live).
# Load the same function dependencies used by the final Figure 4 simulation.
source("01_seeding_functions.r")
source("02_adoption_model.r")
source("03_estimators.r")
source("04_simulate_experiment.r")
source("05_monte_carlo.r")
source("06_load_cai_networks.r")
source("12_optimized_simulation.r")

set.seed(6165)

input_dir <- "results/extension3/clustering_pairs"

output_dir <- "results/extension3/clustering_simulation"

checkpoint_dir <- file.path(output_dir, "checkpoints")

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Representative setting from our original extension
alpha <- -2
beta <- 2
gamma <- 0
rho <- 0.5

n_reps <- 500

# Check that existing simulation functions are loaded
stopifnot(
  exists("prepare_network", mode = "function"),
  exists("simulate_experiment_fast", mode = "function")
)


# ------------------------------------------------------------
# 2. Load the matched network populations
# ------------------------------------------------------------

pools <- readRDS(
  file.path(input_dir, "clustering_network_pools.rds")
)

low_graphs <- pools$low
high_graphs <- pools$high

stopifnot(
  length(low_graphs) == 150,
  length(high_graphs) == 150
)

# Verify each node's degree matches within every pair
for (i in seq_along(low_graphs)) {
  stopifnot(
    identical(
      igraph::degree(low_graphs[[i]]),
      igraph::degree(high_graphs[[i]])
    )
  )
}

cat("Loaded 150 matched network pairs.\n")


# ------------------------------------------------------------
# 3. Prepare networks using our existing function
# ------------------------------------------------------------

cat("Preparing low-clustering networks...\n")

low_prepared <- lapply(
  low_graphs,
  prepare_network
)

cat("Preparing high-clustering networks...\n")

high_prepared <- lapply(
  high_graphs,
  prepare_network
)

# These prepared objects contain the seed probabilities
# and network information used by our fast simulator.


# ------------------------------------------------------------
# 4. Calculate estimators from one simulated experiment
# ------------------------------------------------------------

calculate_estimators <- function(df) {

  required_cols <- c(
    "strategy",
    "outcome",
    "p_A",
    "p_B",
    "p_design",
    "y_A",
    "y_B"
  )

  stopifnot(all(required_cols %in% names(df)))

  stopifnot(
    all(is.finite(df$outcome)),
    all(is.finite(df$p_design)),
    all(df$p_design > 0)
  )

  # These are the importance weights from our mixed design
  w_A <- df$p_A / df$p_design
  w_B <- df$p_B / df$p_design

  # Our simulator records A and B as strategy labels
  strategy <- toupper(trimws(as.character(df$strategy)))

  if (!all(strategy %in% c("A", "B"))) {
    stop(
      "Unexpected strategy labels: ",
      paste(unique(strategy), collapse = ", "),
      ". Check how the original simulator codes A/B."
    )
  }

  if (!any(strategy == "A") || !any(strategy == "B")) {
    stop("Both strategies must appear for DM.")
  }

  # Difference-in-means
  dm <- mean(df$outcome[strategy == "A"]) -
    mean(df$outcome[strategy == "B"])

  # Horvitz-Thompson
  ht <- mean(w_A * df$outcome) -
    mean(w_B * df$outcome)

  # Hajek: normalize by realized weight sums
  hajek <- sum(w_A * df$outcome) / sum(w_A) -
    sum(w_B * df$outcome) / sum(w_B)

  # Simulated benchmark:
  # mean difference in the two generated potential outcomes
  tau_benchmark <- mean(df$y_A - df$y_B)

  tibble(
    DM = dm,
    HT = ht,
    Hajek = hajek,
    tau_benchmark = tau_benchmark
  )
}


# ------------------------------------------------------------
# 5. Run 500 repeated experiments with checkpoints
# ------------------------------------------------------------

for (rep in seq_len(n_reps)) {

  checkpoint_file <- file.path(
    checkpoint_dir,
    sprintf("rep_%04d.rds", rep)
  )

  # Skip finished repetitions if restarting
  if (file.exists(checkpoint_file)) {
    next
  }

  set.seed(616500 + rep)

  # Low-clustering population
  df_low <- simulate_experiment_fast(
    networks = low_prepared,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    rho = rho
  )

  # High-clustering population
  df_high <- simulate_experiment_fast(
    networks = high_prepared,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    rho = rho
  )

  stopifnot(
    nrow(df_low) == 150,
    nrow(df_high) == 150
  )

  # Experiment-level estimator results
  est_low <- calculate_estimators(df_low) %>%
    mutate(
      rep = rep,
      condition = "Low"
    )

  est_high <- calculate_estimators(df_high) %>%
    mutate(
      rep = rep,
      condition = "High"
    )

  # Village-level simulated potential-outcome contrasts
  # Row i corresponds to the same degree sequence in both groups.
  village_results <- tibble(
    rep = rep,
    network_id = seq_len(150),

    tau_low = df_low$y_A - df_low$y_B,

    tau_high = df_high$y_A - df_high$y_B
  ) %>%
    mutate(
      delta_tau = tau_high - tau_low
    )

  # Save each complete repetition
  saveRDS(
    list(
      estimators = bind_rows(est_low, est_high),
      villages = village_results
    ),
    checkpoint_file
  )

  if (rep %% 25 == 0) {
    cat("Completed", rep, "of", n_reps, "repetitions\n")
  }
}

cat("\nMonte Carlo simulations complete.\n")


# ------------------------------------------------------------
# 6. Combine all checkpoints
# ------------------------------------------------------------

checkpoint_files <- file.path(
  checkpoint_dir,
  sprintf("rep_%04d.rds", seq_len(n_reps))
)

stopifnot(all(file.exists(checkpoint_files)))

results <- lapply(checkpoint_files, readRDS)

estimator_results <- bind_rows(
  lapply(results, `[[`, "estimators")
)

village_results <- bind_rows(
  lapply(results, `[[`, "villages")
)

saveRDS(
  estimator_results,
  file.path(output_dir, "estimator_results.rds")
)

saveRDS(
  village_results,
  file.path(output_dir, "village_results.rds")
)


# ============================================================
# ANALYSIS A: TREATMENT EFFECTS
# ============================================================


# ------------------------------------------------------------
# 7. Average treatment effect for each network pair
# ------------------------------------------------------------

paired_effects <- village_results %>%
  group_by(network_id) %>%
  summarise(
    tau_low = mean(tau_low),
    tau_high = mean(tau_high),
    delta_tau = mean(delta_tau),
    .groups = "drop"
  )

# Add the network structural changes
paired_metrics <- read_csv(
  file.path(
    input_dir,
    "clustering_paired_metrics.csv"
  ),
  show_col_types = FALSE
)

paired_effects <- paired_effects %>%
  left_join(
    paired_metrics %>%
      select(
        network_id,
        delta_clustering,
        delta_path_length,
        delta_assortativity
      ),
    by = "network_id"
  )


# ------------------------------------------------------------
# 8. Summarize paired treatment-effect differences
# ------------------------------------------------------------

effect_summary <- paired_effects %>%
  summarise(
    mean_tau_low = mean(tau_low),
    mean_tau_high = mean(tau_high),

    mean_delta_tau = mean(delta_tau),

    mean_absolute_delta = mean(abs(delta_tau)),

    median_delta_tau = median(delta_tau),

    proportion_positive = mean(delta_tau > 0),

    sd_delta_tau = sd(delta_tau)
  )

cat("\nTreatment-effect comparison:\n")
print(effect_summary, width = Inf)


# ------------------------------------------------------------
# 9. Plot paired treatment-effect changes
# ------------------------------------------------------------

p_effects <- ggplot(
  paired_effects,
  aes(x = delta_tau)
) +
  geom_histogram(
    bins = 25,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Effect of Degree-Preserving Topology Changes",
    subtitle = "150 matched high- and low-clustering network pairs",
    x = "Treatment-effect difference (high minus low)",
    y = "Number of network pairs"
  ) +
  theme_minimal(base_size = 12)

print(p_effects)

ggsave(
  file.path(output_dir, "paired_treatment_effects.png"),
  p_effects,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 10. Relate structural changes to treatment-effect changes
# ------------------------------------------------------------

structural_relationships <- paired_effects %>%
  select(
    delta_tau,
    delta_clustering,
    delta_path_length,
    delta_assortativity
  ) %>%
  cor(use = "complete.obs")

cat("\nCorrelations among changes:\n")
print(structural_relationships)

p_structure <- ggplot(
  paired_effects,
  aes(
    x = delta_clustering,
    y = delta_tau
  )
) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    title = "Clustering Changes and Seeding Effectiveness",
    subtitle = "Matched networks with identical degree sequences",
    x = "Clustering difference (high minus low)",
    y = "Treatment-effect difference (high minus low)"
  ) +
  theme_minimal(base_size = 12)

print(p_structure)

ggsave(
  file.path(output_dir, "clustering_vs_treatment_effect.png"),
  p_structure,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ============================================================
# ANALYSIS B: ESTIMATOR EFFICIENCY
# ============================================================


# ------------------------------------------------------------
# 11. Convert estimator results to long format
# ------------------------------------------------------------

estimator_long <- estimator_results %>%
  select(
    rep,
    condition,
    DM,
    HT,
    Hajek
  ) %>%
  pivot_longer(
    cols = c(DM, HT, Hajek),
    names_to = "estimator",
    values_to = "estimate"
  )


# ------------------------------------------------------------
# 12. Compare empirical estimator variances
# ------------------------------------------------------------

variance_summary <- estimator_long %>%
  group_by(condition, estimator) %>%
  summarise(
    mean_estimate = mean(estimate),
    empirical_variance = var(estimate),
    empirical_sd = sd(estimate),
    .groups = "drop"
  )

cat("\nEstimator variance summary:\n")
print(variance_summary, n = Inf)


# ------------------------------------------------------------
# 13. Calculate relative efficiency
# ------------------------------------------------------------

efficiency_summary <- variance_summary %>%
  select(
    condition,
    estimator,
    empirical_variance
  ) %>%
  pivot_wider(
    names_from = estimator,
    values_from = empirical_variance
  ) %>%
  mutate(
    RE_Hajek = DM / Hajek,
    RE_HT = DM / HT
  )

cat("\nRelative efficiency comparison:\n")
print(efficiency_summary, width = Inf)


# ------------------------------------------------------------
# 14. Plot estimator variances
# ------------------------------------------------------------

p_variance <- ggplot(
  variance_summary,
  aes(
    x = estimator,
    y = empirical_variance,
    fill = condition
  )
) +
  geom_col(
    position = "dodge",
    width = 0.65
  ) +
  labs(
    title = "Estimator Variance by Network Topology",
    subtitle = "500 Monte Carlo experiments per condition",
    x = "Estimator",
    y = "Empirical variance",
    fill = "Clustering"
  ) +
  theme_minimal(base_size = 12)

print(p_variance)

ggsave(
  file.path(output_dir, "estimator_variance_comparison.png"),
  p_variance,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 15. Save summaries
# ------------------------------------------------------------

write_csv(
  paired_effects,
  file.path(output_dir, "paired_treatment_effects.csv")
)

write_csv(
  effect_summary,
  file.path(output_dir, "treatment_effect_summary.csv")
)

write_csv(
  variance_summary,
  file.path(output_dir, "estimator_variance_summary.csv")
)

write_csv(
  efficiency_summary,
  file.path(output_dir, "relative_efficiency_summary.csv")
)

cat("\nAll analyses completed.\n")
cat("Results saved to:", output_dir, "\n")






# ------------------------------------------------------------
# 16. Bootstrap uncertainty in relative efficiency difference
# ------------------------------------------------------------

set.seed(6165)

# One row per repetition, with both conditions' estimates
paired_reps <- estimator_results %>%
  select(rep, condition, DM, Hajek) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(DM, Hajek),
    names_sep = "_"
  )

stopifnot(nrow(paired_reps) == 500)

re_value <- function(df) {
  re_low <- var(df$DM_Low) / var(df$Hajek_Low)
  re_high <- var(df$DM_High) / var(df$Hajek_High)

  c(
    RE_low = re_low,
    RE_high = re_high,
    RE_difference = re_high - re_low
  )
}

n_boot <- 2000

boot_results <- replicate(n_boot, {
  ids <- sample(
    seq_len(nrow(paired_reps)),
    replace = TRUE
  )

  re_value(paired_reps[ids, ])
})

# Original estimate and percentile bootstrap intervals
bootstrap_summary <- tibble(
  statistic = rownames(boot_results),
  estimate = unname(re_value(paired_reps)),
  lower_95 = apply(
    boot_results,
    1,
    quantile,
    probs = 0.025
  ),
  upper_95 = apply(
    boot_results,
    1,
    quantile,
    probs = 0.975
  )
)

print(bootstrap_summary)

write_csv(
  bootstrap_summary,
  file.path(
    output_dir,
    "relative_efficiency_bootstrap.csv"
  )
)





# ------------------------------------------------------------
# 17. Bootstrap CI for mean paired treatment-effect difference
# ------------------------------------------------------------

set.seed(6165)

n_boot <- 5000

# One estimated treatment-effect difference per matched pair
paired_differences <- paired_effects$delta_tau

stopifnot(
  length(paired_differences) == 150,
  all(is.finite(paired_differences))
)

# Resample entire matched pairs, keeping their differences intact
bootstrap_means <- replicate(n_boot, {

  sampled_pairs <- sample(
    paired_differences,
    size = length(paired_differences),
    replace = TRUE
  )

  mean(sampled_pairs)
})

# Original point estimate
mean_difference <- mean(paired_differences)

# Percentile bootstrap 95% interval
ci <- quantile(
  bootstrap_means,
  probs = c(0.025, 0.975)
)

bootstrap_effect_summary <- tibble(
  mean_difference = mean_difference,
  lower_95 = unname(ci[1]),
  upper_95 = unname(ci[2])
)

cat("\nBootstrap CI for paired treatment-effect difference:\n")
print(bootstrap_effect_summary, width = Inf)

write_csv(
  bootstrap_effect_summary,
  file.path(
    output_dir,
    "paired_treatment_effect_bootstrap.csv"
  )
)





# Create the comparison plot directly from the calculated treatment effects.
# effect_summary stores adoption proportions, so multiply by 100 for pp.
effect_plot_data <- tibble(
  condition = c("Low clustering", "High clustering"),
  treatment_effect = 100 * c(
    effect_summary$mean_tau_low,
    effect_summary$mean_tau_high
  )
)

p_effect_comparison <- ggplot(
  effect_plot_data,
  aes(x = condition, y = treatment_effect, fill = condition)
) +
  geom_col(width = 0.55) +
  geom_text(
    aes(label = sprintf("%.1f pp", treatment_effect)),
    vjust = -0.5,
    size = 5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "One-Hop Targeting Advantage by Network Topology",
    subtitle = "150 matched network pairs with identical node degrees",
    x = NULL,
    y = "Treatment effect (percentage points)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  file.path(output_dir, "treatment_effect_comparison.png"),
  p_effect_comparison,
  width = 8,
  height = 5,
  dpi = 300
)
