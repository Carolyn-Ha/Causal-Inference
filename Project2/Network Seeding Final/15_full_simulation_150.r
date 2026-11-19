
# ==================================================
# 1. Load packages and simulation functions
# ==================================================

library(tidyverse)
library(igraph)

source("01_seeding_functions.r")
source("02_adoption_model.r")
source("03_estimators.r")
source("04_simulate_experiment.r")
source("05_monte_carlo.r")
source("06_load_cai_networks.r")
source("12_optimized_simulation.r")

# Prepare the Cai networks explicitly, now that file 12 only
# defines functions and no longer creates networks_fast.
networks_fast <- lapply(networks_cai, prepare_network)
names(networks_fast) <- names(networks_cai)
stopifnot(length(networks_fast) >= 150)

# Create output directory
dir.create(
  "results/figure4_full_150",
  recursive = TRUE,
  showWarnings = FALSE
)




# ==================================================
# 2. Full Figure 4 parameter grid
# ==================================================

params_full_figure4 <- expand.grid(
  alpha = c(-3, -2, -1, 0),
  beta = seq(0, 8, by = 2),
  gamma = c(0, 0.1),
  n_villages = 150,
  rho = c(0, 0.5)
) %>%
  as_tibble() %>%
  mutate(setting_id = row_number())

# Save the parameter grid
write_csv(
  params_full_figure4,
  "results/figure4_full_150/parameter_grid.csv"
)

nrow(params_full_figure4)




# ==================================================
# 3. Checkpointed Monte Carlo runner
# ==================================================

run_checkpointed_figure4 <- function(
    params,
    networks,
    R = 5000,
    chunk_size = 100,
    output_dir = "results/figure4_full_150"
) {

  if (R %% chunk_size != 0) {
    stop("R must be divisible by chunk_size.")
  }

  n_chunks <- R / chunk_size

  for (i in seq_len(nrow(params))) {

    p <- params[i, ]

    setting_dir <- file.path(
      output_dir,
      sprintf("setting_%03d", p$setting_id)
    )

    dir.create(
      setting_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    cat(
      "\nSetting", i, "of", nrow(params),
      "| alpha =", p$alpha,
      "| beta =", p$beta,
      "| gamma =", p$gamma,
      "| rho =", p$rho, "\n"
    )

    for (chunk in seq_len(n_chunks)) {

      checkpoint_file <- file.path(
        setting_dir,
        sprintf("chunk_%03d.rds", chunk)
      )

      # Skip completed checkpoints
      if (file.exists(checkpoint_file)) {
        next
      }

      # Reproducible seed for this setting and chunk
      set.seed(
        6165000 + p$setting_id * 1000 + chunk
      )

      # Simulate this chunk
      chunk_results <- run_monte_carlo_fast(
        networks = networks,
        alpha = p$alpha,
        beta = p$beta,
        gamma = p$gamma,
        R = chunk_size,
        n_villages = p$n_villages,
        rho = p$rho
      )

      # Adjust repetition numbers so they are unique
      # across chunks of the same setting
      chunk_results <- chunk_results %>%
        mutate(
          rep = rep + (chunk - 1) * chunk_size,
          setting_id = p$setting_id
        )

      # Save through a temporary file first
      temp_file <- tempfile(
        pattern = "checkpoint_",
        tmpdir = setting_dir,
        fileext = ".rds"
      )

      saveRDS(chunk_results, temp_file)

      if (!file.rename(temp_file, checkpoint_file)) {
        stop("Could not save: ", checkpoint_file)
      }

      cat(
        "  Saved chunk", chunk, "of", n_chunks,
        "\n"
      )
    }
  }

  cat("\nAll simulation checkpoints complete!\n")
}



# ==================================================
# 4. Run full simulation
# ==================================================

run_checkpointed_figure4(
  params = params_full_figure4,
  networks = networks_fast,
  R = 5000,
  chunk_size = 100,
  output_dir = "results/figure4_full_150"
)