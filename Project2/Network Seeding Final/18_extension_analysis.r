
# ==================================================
# 18_extension_analysis.r
# Extension 3: Main treatment-effect and efficiency analysis
# Step 1. Load and summarize simulation results
# ==================================================

library(tidyverse)

# Prefer saved combined results; reconstruct from checkpoints if needed.
combined_file <- "results/extension3/extension_results_full.rds"

if (file.exists(combined_file)) {
  extension_results <- readRDS(combined_file)
} else {
  extension_files <- list.files(
    "results/extension3/main_simulation",
    pattern = "^chunk_[0-9]+\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  # 60 settings x 20 chunks of 100 repetitions.
  stopifnot(length(extension_files) == 60L * 20L)
  extension_results <- purrr::map_dfr(extension_files, readRDS)
  saveRDS(extension_results, combined_file)
}

# Each of the 60 settings has 2,000 distinct repetitions per estimator.
repetition_checks <- extension_results %>%
  group_by(attachment, alpha, beta, gamma, estimator) %>%
  summarise(
    n_reps = n(),
    distinct_reps = n_distinct(rep),
    .groups = "drop"
  )

stopifnot(
  nrow(repetition_checks) == 60L * 3L,
  all(repetition_checks$n_reps == 2000L),
  all(repetition_checks$distinct_reps == 2000L)
)

# Summarize estimator performance
extension_summary <- extension_results %>%
  group_by(
    attachment, alpha, beta, gamma, estimator
  ) %>%
  summarise(
    n_reps = n(),
    mean_estimate = mean(estimate),
    empirical_variance = var(estimate),
    empirical_SD = sd(estimate),
    mean_SE = mean(SE),
    rejection_rate = mean(reject_H0),
    .groups = "drop"
  )

# Save summary
write_csv(
  extension_summary,
  "results/extension3/extension_summary.csv"
)

print(extension_summary, n = 15, width = Inf)



# ==================================================
# Step 2. Treatment effect across attachment levels
# ==================================================

# The simulation records both potential outcomes within each repetition.
# Remove repeated estimator rows before averaging the benchmark effect.
treatment_effect_summary <- extension_results %>%
  distinct(
    attachment, alpha, beta, gamma, rep, tau_benchmark
  ) %>%
  group_by(attachment, alpha, beta, gamma) %>%
  summarise(
    n_reps = n(),
    tau_benchmark = mean(tau_benchmark),
    .groups = "drop"
  ) %>%
  filter(alpha == -2, beta == 2, gamma == 0) %>%
  arrange(attachment)

stopifnot(all(treatment_effect_summary$n_reps == 2000L))

print(treatment_effect_summary, width = Inf)

# Plot the simulated benchmark (not the mean DM estimate).
ggplot(
  treatment_effect_summary,
  aes(x = attachment, y = tau_benchmark)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2)
  ) +
  labs(
    x = "Preferential attachment parameter",
    y = "Simulated benchmark treatment effect",
    title = "Effect of network heterogeneity on one-hop targeting",
    subtitle = expression(
      alpha == -2 ~ "," ~ beta == 2 ~ "," ~ gamma == 0
    )
  ) +
  theme_bw(base_size = 12)

# ==================================================
# Step 3. Relative efficiency across attachment levels
# ==================================================

efficiency_summary <- extension_summary %>%
  filter(
    alpha == -2,
    beta == 2,
    gamma == 0
  ) %>%
  select(
    attachment,
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

print(efficiency_summary, width = Inf)

# Plot efficiency ratios
efficiency_long <- efficiency_summary %>%
  select(attachment, RE_Hajek, RE_HT) %>%
  pivot_longer(
    cols = c(RE_Hajek, RE_HT),
    names_to = "estimator",
    values_to = "relative_efficiency"
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      RE_Hajek = "Hajek",
      RE_HT = "Horvitz-Thompson"
    )
  )

ggplot(
  efficiency_long,
  aes(
    x = attachment,
    y = relative_efficiency,
    color = estimator
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 1.5, 2)
  ) +
  labs(
    x = "Preferential attachment parameter",
    y = "Relative efficiency (DM variance / estimator variance)",
    color = "Estimator",
    title = "Estimator efficiency across network structures"
  ) +
  theme_bw(base_size = 12)




# ==================================================
# Step 4. Relative efficiency across all adoption settings
# ==================================================

efficiency_all <- extension_summary %>%
  select(
    attachment, alpha, beta, gamma,
    estimator, empirical_variance
  ) %>%
  pivot_wider(
    names_from = estimator,
    values_from = empirical_variance
  ) %>%
  mutate(
    RE_Hajek = DM / Hajek,
    RE_HT = DM / HT
  )

# Plot Hájek relative efficiency across all 12 settings
ggplot(
  efficiency_all,
  aes(
    x = attachment,
    y = RE_Hajek
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(
    gamma + alpha ~ beta,
    labeller = label_both
  ) +
  scale_x_continuous(
    breaks = c(0, 1, 2)
  ) +
  labs(
    x = "Preferential attachment parameter",
    y = "Relative efficiency (DM variance / Hájek variance)",
    title = "Hájek efficiency across network structures",
    subtitle = "Comparison across adoption-model parameters"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

