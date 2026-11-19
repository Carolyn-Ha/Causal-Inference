
# ==================================================
# 17_extension_main_simulation.r
# Extension 3: Main preferential-attachment Monte Carlo simulation
# ==================================================

library(tidyverse)
library(igraph)

# Load existing simulation functions
source("05_monte_carlo.r")
source("12_optimized_simulation.r")

# Load our prepared synthetic networks
network_pools_fast <- readRDS(
  "results/extension3/network_pools_fast.rds"
)

# Define experimental parameter grid
extension_params <- expand.grid(
  attachment = c(0, 0.5, 1, 1.5, 2),
  alpha = c(-2, -1),
  beta = c(0, 2, 4),
  gamma = c(0, 0.1),
  n_villages = 150,
  rho = 0.5
) %>%
  as_tibble() %>%
  mutate(setting_id = row_number())

# Create a separate directory for extension results
dir.create(
  "results/extension3/main_simulation",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save parameter grid
write_csv(
  extension_params,
  "results/extension3/main_simulation/parameter_grid.csv"
)

# Verify
nrow(extension_params)

print(extension_params, n = 10)




# ==================================================
# Step 2. Checkpointed extension simulation
# ==================================================

run_extension_checkpointed <- function(
    params,
    network_pools,
    R = 2000,
    chunk_size = 100,
    output_dir = "results/extension3/main_simulation"
) {

  if (R %% chunk_size != 0) {
    stop("R must be divisible by chunk_size.")
  }

  n_chunks <- R / chunk_size

  for (i in seq_len(nrow(params))) {

    p <- params[i, ]

    # Select the appropriate network pool
    pool_name <- paste0("attachment_", p$attachment)

    networks <- network_pools[[pool_name]]

    if (is.null(networks)) {
      stop("Network pool not found: ", pool_name)
    }

    # Separate folder for each parameter setting
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
      "| attachment =", p$attachment,
      "| alpha =", p$alpha,
      "| beta =", p$beta,
      "| gamma =", p$gamma, "\n"
    )

    for (chunk in seq_len(n_chunks)) {

      checkpoint_file <- file.path(
        setting_dir,
        sprintf("chunk_%03d.rds", chunk)
      )

      # Skip completed chunks
      if (file.exists(checkpoint_file)) {
        next
      }

      # Reproducible seed
      set.seed(
        6170000 + p$setting_id * 1000 + chunk
      )

      # Run Monte Carlo simulations
      chunk_results <- run_monte_carlo_fast(
        networks = networks,
        alpha = p$alpha,
        beta = p$beta,
        gamma = p$gamma,
        R = chunk_size,
        n_villages = p$n_villages,
        rho = p$rho
      )

      # Add metadata and unique repetition numbers
      chunk_results <- chunk_results %>%
        mutate(
          attachment = p$attachment,
          rep = rep + (chunk - 1) * chunk_size,
          setting_id = p$setting_id
        )

      # Save checkpoint safely
      temp_file <- tempfile(
        pattern = "checkpoint_",
        tmpdir = setting_dir,
        fileext = ".rds"
      )

      saveRDS(chunk_results, temp_file)

      if (!file.rename(temp_file, checkpoint_file)) {
        stop("Could not save checkpoint: ", checkpoint_file)
      }

      cat(
        "  Saved chunk", chunk, "of", n_chunks,
        "\n"
      )
    }
  }

  cat("\nAll extension checkpoints complete!\n")
}




# ==================================================
# Step 3. Run full extension simulation
# ==================================================

run_extension_checkpointed(
  params = extension_params,
  network_pools = network_pools_fast,
  R = 2000,
  chunk_size = 100,
  output_dir = "results/extension3/main_simulation"
)

