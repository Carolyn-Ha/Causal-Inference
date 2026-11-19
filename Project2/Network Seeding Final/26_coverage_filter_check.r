# ==================================================
# 26_coverage_filter_check.r
# Coverage sensitivity analysis for the final Figure 4 replication
# ==================================================

# ==================================================
# 1. Load packages and saved simulation results
# ==================================================

library(tidyverse)
library(patchwork)

results <- readRDS(
  "results/figure4_results_150_5000.rds"
)

output_dir <- "results/figure4_coverage"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Candidate threshold from the authors' commented-out code.
# This is a sensitivity analysis, not a verified
# reproduction of their published exclusion criterion.
coverage_cutoff <- 0.85

# Inspect the available results
print(table(results$estimator, results$rho))


# ==================================================
# 2. Approximate the true treatment effect
# ==================================================

# Each simulated experiment's benchmark is repeated
# across estimator rows. Count it only once.

tau_truth <- results %>%
  distinct(
    alpha, beta, gamma, n_villages,
    rho, rep, tau_benchmark
  ) %>%
  group_by(alpha, beta, gamma, n_villages) %>%
  summarise(
    tau_true = mean(tau_benchmark, na.rm = TRUE),
    .groups = "drop"
  )


# ==================================================
# 3. Calculate 90% confidence-interval coverage
# ==================================================

results_with_coverage <- results %>%
  left_join(
    tau_truth,
    by = c("alpha", "beta", "gamma", "n_villages")
  ) %>%
  mutate(
    valid = is.finite(CI_lower) &
            is.finite(CI_upper) &
            is.finite(tau_true),

    covered = case_when(
      valid ~ CI_lower <= tau_true &
              CI_upper >= tau_true,
      .default = NA
    )
  )

coverage_summary <- results_with_coverage %>%
  group_by(
    alpha, beta, gamma, n_villages,
    rho, estimator
  ) %>%
  summarise(
    n_reps = n(),
    n_valid = sum(valid),
    coverage = if (any(valid)) {
      mean(covered[valid])
    } else {
      NA_real_
    },
    power = if (all(!is.na(reject_H0))) {
      mean(reject_H0)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(
    passes_coverage = !is.na(coverage) &
                      coverage >= coverage_cutoff
  )

write_csv(
  coverage_summary,
  file.path(output_dir, "coverage_summary.csv")
)

# Check whether any estimator-setting combination has invalid CIs.
# If every interval is valid, this prints zero rows.
cat("\nSettings with invalid confidence intervals:\n")
coverage_summary %>%
  filter(n_valid < n_reps) %>%
  print(n = Inf)


# ==================================================
# 4. Examine the coverage results
# ==================================================

coverage_diagnostics <- coverage_summary %>%
  group_by(estimator) %>%
  summarise(
    total_settings = n(),
    valid_settings = sum(!is.na(coverage)),
    below_85 = sum(
      !is.na(coverage) & coverage < coverage_cutoff
    ),
    mean_coverage = mean(coverage, na.rm = TRUE),
    min_coverage = min(coverage, na.rm = TRUE),
    .groups = "drop"
  )

print(coverage_diagnostics, n = Inf)

# Show every combination below the threshold or with undefined coverage
low_coverage <- coverage_summary %>%
  filter(!passes_coverage) %>%
  arrange(coverage)

print(low_coverage, n = Inf)

write_csv(
  low_coverage,
  file.path(output_dir, "low_coverage_settings.csv")
)


# ==================================================
# 5. Prepare data for Figure 4
# ==================================================

estimator_levels <- c(
  "Difference-in-means",
  "Horvitz-Thompson",
  "Hajek",
  "HT (off-policy)",
  "Hajek (off-policy)"
)

estimator_colors <- c(
  "Difference-in-means" = "#E69F00",
  "Horvitz-Thompson" = "#56B4E9",
  "Hajek" = "#009E73",
  "HT (off-policy)" = "#1F78B4",
  "Hajek (off-policy)" = "#007F5F"
)

estimator_shapes <- c(
  "Difference-in-means" = 0,
  "Horvitz-Thompson" = 0,
  "Hajek" = 0,
  "HT (off-policy)" = 4,
  "Hajek (off-policy)" = 4
)

plot_data <- coverage_summary %>%
  left_join(
    tau_truth,
    by = c("alpha", "beta", "gamma", "n_villages")
  ) %>%
  mutate(
    estimator = recode(
      estimator,
      "DM" = "Difference-in-means",
      "HT" = "Horvitz-Thompson",
      "Hajek" = "Hajek",
      "HT_offpolicy" = "HT (off-policy)",
      "Hajek_offpolicy" = "Hajek (off-policy)"
    ),
    estimator = factor(
      estimator,
      levels = estimator_levels
    )
  )

# Exclude beta=gamma=0 as in our original figure.
left_data_all <- plot_data %>%
  filter(beta > 0 | gamma > 0)

# Retain *all* settings but set failed coverage values to NA.
# This makes excluded points explicit in the plotting data.
left_data_filtered <- left_data_all %>%
  mutate(
    power_filtered = if_else(
      passes_coverage,
      power,
      NA_real_
    )
  )

cat(
  "Left-panel points before filtering:",
  nrow(left_data_all), "\n"
)

cat(
  "Left-panel points after filtering:",
  sum(!is.na(left_data_filtered$power_filtered)), "\n"
)

cat(
  "Points removed:",
  sum(is.na(left_data_filtered$power_filtered)), "\n"
)

# Important: geom_line() can bridge across missing values.
# To ensure coverage exclusions create visible gaps, build
# explicit segments ONLY between adjacent retained beta settings.
# The original beta grid is 0, 2, 4, 6, 8.
left_line_segments <- left_data_filtered %>%
  group_by(alpha, gamma, n_villages, rho, estimator) %>%
  arrange(beta, .by_group = TRUE) %>%
  mutate(
    beta_previous = lag(beta),
    power_previous = lag(power_filtered)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(power_filtered),
    !is.na(power_previous),
    beta - beta_previous == 2
  )


# ==================================================
# 6. Create coverage-filtered LEFT panels
# ==================================================

figure4_left_filtered <- ggplot(
  left_data_filtered,
  aes(
    x = beta,
    y = power_filtered,
    color = estimator,
    shape = estimator,
    group = estimator
  )
) +
  geom_hline(
    yintercept = 0.10,
    linetype = "dashed",
    color = "black",
    linewidth = 0.4
  ) +
  geom_segment(
    data = left_line_segments,
    aes(
      x = beta_previous,
      xend = beta,
      y = power_previous,
      yend = power_filtered,
      color = estimator
    ),
    inherit.aes = FALSE,
    linewidth = 0.75
  ) +
  geom_point(size = 2, stroke = 0.9, na.rm = TRUE) +
  facet_grid(gamma ~ alpha) +
  scale_color_manual(
    values = estimator_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = estimator_shapes,
    drop = FALSE
  ) +
  scale_x_continuous(breaks = seq(0, 8, by = 2)) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    x = expression("Spillover effect " * beta),
    y = "Power"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey85"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.spacing = grid::unit(0.15, "cm")
  )


# ==================================================
# 7. Create upper-right panel
#    Keep unfiltered, as in our previous figure
# ==================================================

figure4_top_right <- plot_data %>%
  filter(
    estimator %in% c(
      "Difference-in-means",
      "Horvitz-Thompson",
      "Hajek"
    ),
    beta > 0 | gamma > 0,
    tau_true > 0
  ) %>%
  ggplot(aes(
    x = tau_true,
    y = power,
    color = estimator,
    shape = estimator
  )) +
  geom_hline(
    yintercept = 0.10,
    linetype = "dotted"
  ) +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_color_manual(values = estimator_colors) +
  scale_shape_manual(values = c(
    "Difference-in-means" = 16,
    "Horvitz-Thompson" = 15,
    "Hajek" = 17
  )) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    x = expression("Treatment effect " * tau),
    y = "Power"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")


# ==================================================
# 8. Create lower-right CDF
#    Keep unfiltered, as in our previous figure
# ==================================================

dm_power <- plot_data %>%
  filter(estimator == "Difference-in-means") %>%
  select(
    alpha, beta, gamma, n_villages,
    power_dm = power
  )

power_differences <- plot_data %>%
  filter(
    estimator %in% c("Horvitz-Thompson", "Hajek"),
    beta > 0 | gamma > 0
  ) %>%
  left_join(
    dm_power,
    by = c("alpha", "beta", "gamma", "n_villages")
  ) %>%
  mutate(
    power_difference = power - power_dm
  )

figure4_bottom_right <- ggplot(
  power_differences,
  aes(x = power_difference, color = estimator)
) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  scale_color_manual(values = estimator_colors) +
  labs(
    x = "Difference in power",
    y = "CDF"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")


# ==================================================
# 9. Combine and save new Figure 4
# ==================================================

figure4_coverage_filtered <- (
  figure4_left_filtered |
  (figure4_top_right / figure4_bottom_right)
) +
  plot_layout(widths = c(2.3, 1))

print(figure4_coverage_filtered)

ggsave(
  file.path(
    output_dir,
    "figure4_150_coverage_filtered.png"
  ),
  figure4_coverage_filtered,
  width = 16,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "figure4_150_coverage_filtered.pdf"
  ),
  figure4_coverage_filtered,
  width = 16,
  height = 8
)

cat("\nCoverage analysis completed.\n")

