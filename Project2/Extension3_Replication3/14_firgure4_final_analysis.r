# ============================================================
# 14_figure4_final_analysis.r
# Final Figure 4: 150-village, 5000-repetition checkpoint analysis
# Requires completed checkpoints from 15_full_simulation_150.r
# Does NOT rerun the Monte Carlo simulation.
# ============================================================

library(tidyverse)
library(patchwork)

# Estimator labels and plotting conventions
estimator_levels <- c(
  "Difference-in-means", "Horvitz-Thompson", "Hajek",
  "HT (off-policy)", "Hajek (off-policy)"
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

dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ==================================================
# 150-village simulation: Load results
# ==================================================

combined_file <- "results/figure4_results_150_5000.rds"

if (file.exists(combined_file)) {

  # Use previously combined simulation results
  figure4_results_150 <- readRDS(combined_file)

} else {

  # Reconstruct results from simulation checkpoints
  checkpoint_files_150 <- list.files(
    "results/figure4_full_150",
    pattern = "^chunk_[0-9]+\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  stopifnot(length(checkpoint_files_150) == 4000)

  figure4_results_150 <- purrr::map_dfr(
    checkpoint_files_150,
    readRDS
  )

  # Save combined results for future use
  saveRDS(figure4_results_150, combined_file)
}

# Validate the loaded results
stopifnot(
  nrow(figure4_results_150) == 1000000,
  all(figure4_results_150$n_villages == 150)
)

expected_estimator_counts <- figure4_results_150 %>%
  count(setting_id, estimator, name = "n_reps")

stopifnot(
  nrow(expected_estimator_counts) == 200,
  all(expected_estimator_counts$n_reps == 5000)
)





# ==================================================
# 150-village simulation: Summarize results
# ==================================================

# Calculate rejection probabilities for each estimator
figure4_power_150 <- figure4_results_150 %>%
  group_by(alpha, beta, gamma, n_villages, estimator) %>%
  summarise(
    n_reps = n(),
    valid_tests = sum(!is.na(reject_H0)),
    power = if (all(!is.na(reject_H0))) {
      mean(reject_H0)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# Calculate average treatment effect for each setting
figure4_tau_150 <- figure4_results_150 %>%
  distinct(
    alpha, beta, gamma, n_villages,
    rho, rep, tau_benchmark
  ) %>%
  group_by(alpha, beta, gamma, n_villages) %>%
  summarise(
    tau = mean(tau_benchmark, na.rm = TRUE),
    .groups = "drop"
  )

# Combine power and treatment effects
figure4_summary_150 <- figure4_power_150 %>%
  left_join(
    figure4_tau_150,
    by = c("alpha", "beta", "gamma", "n_villages")
  )

# Save summarized results
saveRDS(
  figure4_summary_150,
  "results/figure4_summary_150_5000.rds"
)

write_csv(
  figure4_summary_150,
  "results/figure4_summary_150_5000.csv"
)

# Verify summary
nrow(figure4_summary_150)

print(head(figure4_summary_150), width = Inf)



# ==================================================
# 150-village simulation: Prepare plotting labels
# ==================================================

figure4_summary_150 <- figure4_summary_150 %>%
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

# Check labels
table(figure4_summary_150$estimator, useNA = "always")





# ==================================================
# 150-village Figure 4: Left panels
# ==================================================

figure4_left_150 <- ggplot(
  figure4_summary_150 %>%
    filter(beta > 0 | gamma > 0),
  aes(
    x = beta,
    y = power,
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
  geom_line(linewidth = 0.75) +
  geom_point(size = 2, stroke = 0.9) +
  facet_grid(
    gamma ~ alpha,
    labeller = labeller(
      alpha = function(x) paste0("alpha = ", x),
      gamma = function(x) paste0("gamma = ", x)
    )
  ) +
  scale_color_manual(
    values = estimator_colors,
    name = "Estimator"
  ) +
  scale_shape_manual(
    values = estimator_shapes,
    name = "Estimator"
  ) +
  scale_x_continuous(
    breaks = seq(0, 8, by = 2)
  ) +
  coord_cartesian(ylim = c(0, 1.02)) +
  labs(
    x = expression("Spillover effect " * beta),
    y = "Rejection probability"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.width = grid::unit(1.1, "cm"),
    panel.spacing = grid::unit(0.15, "cm")
  )

print(figure4_left_150)



# ==================================================
# 150-village Figure 4: Upper-right panel
# ==================================================

figure4_top_right_150 <- figure4_summary_150 %>%
  filter(
    estimator %in% c(
      "Difference-in-means",
      "Horvitz-Thompson",
      "Hajek"
    ),
    beta > 0 | gamma > 0,
    tau > 0
  ) %>%
  ggplot(
    aes(
      x = tau,
      y = power,
      color = estimator,
      shape = estimator
    )
  ) +
  geom_hline(
    yintercept = 0.10,
    linetype = "dotted",
    linewidth = 0.5
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
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(figure4_top_right_150)




# ==================================================
# 150-village Figure 4: Lower-right CDF
# ==================================================

# Get difference-in-means power for each setting
dm_power_150 <- figure4_summary_150 %>%
  filter(estimator == "Difference-in-means") %>%
  select(
    alpha, beta, gamma, n_villages,
    power_dm = power
  )

# Calculate differences in power relative to DM
power_differences_150 <- figure4_summary_150 %>%
  filter(
    estimator %in% c(
      "Horvitz-Thompson",
      "Hajek"
    ),
    beta > 0 | gamma > 0
  ) %>%
  left_join(
    dm_power_150,
    by = c("alpha", "beta", "gamma", "n_villages")
  ) %>%
  mutate(
    power_difference = power - power_dm
  )

# Plot empirical CDF
figure4_bottom_right_150 <- ggplot(
  power_differences_150,
  aes(
    x = power_difference,
    color = estimator
  )
) +
  stat_ecdf(
    geom = "step",
    linewidth = 0.9
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = estimator_colors
  ) +
  labs(
    x = "Difference in power",
    y = "CDF"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(figure4_bottom_right_150)




# ==================================================
# Final 150-village Figure 4
# ==================================================

figure4_complete_150 <- (
  figure4_left_150 |
    (
      figure4_top_right_150 /
      figure4_bottom_right_150
    )
) +
  plot_layout(widths = c(2.3, 1))

print(figure4_complete_150)

ggsave(
  "figures/figure4_final_150.png",
  plot = figure4_complete_150,
  width = 16,
  height = 8,
  dpi = 300
)

ggsave(
  "figures/figure4_final_150.pdf",
  plot = figure4_complete_150,
  width = 16,
  height = 8
)





