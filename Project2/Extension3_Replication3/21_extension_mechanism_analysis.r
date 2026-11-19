
# ==================================================
# Extension 3: Mechanism analysis
# Step 1. Importance-weight variability
# ==================================================

library(tidyverse)

# The fast simulator is defined here; it no longer prepares Cai data
# or runs simulations automatically when sourced.
source("12_optimized_simulation.r")

# Load original prepared network pools
network_pools_fast <- readRDS(
  "results/extension3/network_pools_fast.rds"
)

attachment_levels <- c(0, 0.5, 1, 1.5, 2)
stopifnot(
  length(network_pools_fast) == length(attachment_levels),
  all(lengths(network_pools_fast) == 150L)
)

set.seed(6172)

weight_diagnostics <- purrr::map2_dfr(
  network_pools_fast,
  attachment_levels,
  function(pool, a) {

    purrr::map_dfr(seq_len(100), function(r) {

      df <- simulate_experiment_fast(
        networks = pool,
        alpha = -2,
        beta = 2,
        gamma = 0,
        rho = 0.5
      )

      w_A <- df$p_A / df$p_design
      w_B <- df$p_B / df$p_design

      tibble(
        attachment = a,
        rep = r,

        # Variation among village weights
        SD_w_A = sd(w_A),
        SD_w_B = sd(w_B),

        # Fluctuation of total weights
        mean_w_A = mean(w_A),
        mean_w_B = mean(w_B),

        # Effective sample sizes
        ESS_A = sum(w_A)^2 / sum(w_A^2),
        ESS_B = sum(w_B)^2 / sum(w_B^2)
      )
    })
  }
)

# Summarize across 100 experiments
weight_summary <- weight_diagnostics %>%
  group_by(attachment) %>%
  summarise(
    avg_SD_w_A = mean(SD_w_A),
    avg_SD_w_B = mean(SD_w_B),

    SD_total_w_A = sd(mean_w_A),
    SD_total_w_B = sd(mean_w_B),

    avg_ESS_A = mean(ESS_A),
    avg_ESS_B = mean(ESS_B),

    .groups = "drop"
  )

print(weight_summary, width = Inf)




# ==================================================
# Step 2. Village-level outcome variability
# ==================================================

set.seed(6173)

outcome_diagnostics <- purrr::map2_dfr(
  network_pools_fast,
  attachment_levels,
  function(pool, a) {

    purrr::map_dfr(seq_len(100), function(r) {

      df <- simulate_experiment_fast(
        networks = pool,
        alpha = -2,
        beta = 2,
        gamma = 0,
        rho = 0.5
      )

      tibble(
        attachment = a,
        rep = r,

        mean_y_A = mean(df$y_A),
        mean_y_B = mean(df$y_B),

        SD_y_A = sd(df$y_A),
        SD_y_B = sd(df$y_B),

        cor_weight_A_y = cor(
          df$p_A / df$p_design,
          df$outcome
        )
      )
    })
  }
)

outcome_summary <- outcome_diagnostics %>%
  group_by(attachment) %>%
  summarise(
    avg_y_A = mean(mean_y_A),
    avg_y_B = mean(mean_y_B),
    avg_SD_y_A = mean(SD_y_A),
    avg_SD_y_B = mean(SD_y_B),
    avg_cor_weight_A_y = mean(cor_weight_A_y),
    .groups = "drop"
  )

print(outcome_summary, width = Inf)



# ==================================================
# Step 3. Examine normalization and weight totals
# ==================================================

set.seed(6174)

normalization_diagnostics <- purrr::map2_dfr(
  network_pools_fast,
  attachment_levels,
  function(pool, a) {

    purrr::map_dfr(seq_len(300), function(r) {

      df <- simulate_experiment_fast(
        networks = pool,
        alpha = -2,
        beta = 2,
        gamma = 0,
        rho = 0.5
      )

      w_A <- df$p_A / df$p_design
      w_B <- df$p_B / df$p_design

      U_A <- mean(w_A * df$outcome)
      U_B <- mean(w_B * df$outcome)

      W_A <- mean(w_A)
      W_B <- mean(w_B)

      tibble(
        attachment = a,
        rep = r,
        HT = U_A - U_B,
        Hajek = U_A / W_A - U_B / W_B,
        W_A = W_A,
        W_B = W_B,
        U_A = U_A,
        U_B = U_B,

        mu_A_HT = U_A,
        mu_B_HT = U_B,

        mu_A_Hajek = U_A / W_A,
        mu_B_Hajek = U_B / W_B

      )
    })
  }
)

# Summarize variation across experiments
normalization_summary <- normalization_diagnostics %>%
  group_by(attachment) %>%
  summarise(
    var_HT = var(HT),
    var_Hajek = var(Hajek),
    variance_ratio_HT_Hajek = var(HT) / var(Hajek),
    SD_W_A = sd(W_A),
    cor_UA_WA = cor(U_A, W_A),
    .groups = "drop"
  )

print(normalization_summary, width = Inf)

# Save the three diagnostic summaries for the report and audit.
output_dir <- "results/extension3/mechanism_analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(weight_summary, file.path(output_dir, "weight_summary.csv"))
write_csv(outcome_summary, file.path(output_dir, "outcome_summary.csv"))
write_csv(normalization_summary, file.path(output_dir, "normalization_summary.csv"))
cat("Mechanism summaries saved to:", output_dir, "\n")

