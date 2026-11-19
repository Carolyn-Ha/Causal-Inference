
# ==================================================
# 19_extension_robustness.r
# Extension 3: Independent-network-pool robustness
# Step 1. Generate additional network populations
# ==================================================

library(tidyverse)
library(igraph)

# File 05 provides summarize_experiment(); File 12 provides
# prepare_network() and run_monte_carlo_fast().
source("05_monte_carlo.r")
source("12_optimized_simulation.r")

set.seed(6171)

attachment_levels <- c(0, 0.5, 1, 1.5, 2)

n_nodes <- 50
m_edges <- 2
n_villages <- 150

# Generate four new independent pools per level
new_pools <- lapply(
  2:5,
  function(pool_id) {

    pools_at_level <- lapply(
      attachment_levels,
      function(a) {

        lapply(
          seq_len(n_villages),
          function(i) {
            sample_pa(
              n = n_nodes,
              power = a,
              m = m_edges,
              directed = FALSE,
              algorithm = "psumtree"
            )
          }
        )
      }
    )

    names(pools_at_level) <- paste0(
      "attachment_", attachment_levels
    )

    pools_at_level
  }
)

names(new_pools) <- paste0("pool_", 2:5)

# Check number of generated networks
pool_counts <- purrr::imap_dfr(
  new_pools,
  function(pool_group, pool_name) {

    tibble(
      pool = pool_name,
      attachment = attachment_levels,
      n_networks = sapply(pool_group, length)
    )
  }
)

print(pool_counts, n = Inf)



# ==================================================
# Step 2. Prepare and save new network pools
# ==================================================

# Convert igraph networks into optimized format
new_pools_fast <- lapply(
  new_pools,
  function(pool_group) {
    lapply(
      pool_group,
      function(pool) {
        lapply(pool, prepare_network)
      }
    )
  }
)

# Save the prepared pools
saveRDS(
  new_pools_fast,
  "results/extension3/robustness_network_pools.rds"
)

robustness_checks <- purrr::imap_dfr(
  new_pools_fast,
  function(pool_group, pool_name) {

    purrr::imap_dfr(
      pool_group,
      function(networks, attachment_name) {

        tibble(
          pool = pool_name,
          attachment = attachment_name,
          n_networks = length(networks),
          min_prob_sum = min(vapply(
            networks,
            function(net) sum(net$q),
            numeric(1)
          )),
          max_prob_sum = max(vapply(
            networks,
            function(net) sum(net$q),
            numeric(1)
          )),
          min_distinct_prob = min(vapply(
            networks,
            function(net) net$pi_distinct,
            numeric(1)
          ))
        )
      }
    )
  }
)

print(robustness_checks, n = Inf)





# ==================================================
# Step 3. Run robustness simulations
# ==================================================

dir.create(
  "results/extension3/robustness_simulation",
  recursive = TRUE,
  showWarnings = FALSE
)

robustness_params <- expand.grid(
  pool = names(new_pools_fast),
  attachment = attachment_levels
) %>%
  as_tibble() %>%
  mutate(setting_id = row_number())

run_robustness <- function(
    params,
    network_pools,
    R = 500,
    chunk_size = 100
) {

  stopifnot(R %% chunk_size == 0)

  for (i in seq_len(nrow(params))) {

    p <- params[i, ]

    pool_name <- paste0("attachment_", p$attachment)

    networks <- network_pools[[p$pool]][[pool_name]]

    if (is.null(networks)) {
      stop("Missing network pool: ", p$pool, " / ", pool_name)
    }

    setting_dir <- file.path(
      "results/extension3/robustness_simulation",
      sprintf("setting_%03d", p$setting_id)
    )

    dir.create(setting_dir, recursive = TRUE, showWarnings = FALSE)

    cat(
      "\nSetting", i, "of", nrow(params),
      "| pool =", p$pool,
      "| attachment =", p$attachment, "\n"
    )

    for (chunk in seq_len(R / chunk_size)) {

      checkpoint_file <- file.path(
        setting_dir,
        sprintf("chunk_%03d.rds", chunk)
      )

      if (file.exists(checkpoint_file)) {
        next
      }

      set.seed(6190000 + p$setting_id * 1000 + chunk)

      chunk_results <- run_monte_carlo_fast(
        networks = networks,
        alpha = -2,
        beta = 2,
        gamma = 0,
        R = chunk_size,
        n_villages = 150,
        rho = 0.5
      ) %>%
        mutate(
          pool = p$pool,
          attachment = p$attachment,
          setting_id = p$setting_id,
          rep = rep + (chunk - 1) * chunk_size
        )

      temp_file <- tempfile(
        pattern = "checkpoint_",
        tmpdir = setting_dir,
        fileext = ".rds"
      )

      saveRDS(chunk_results, temp_file)

      if (!file.rename(temp_file, checkpoint_file)) {
        stop("Could not save checkpoint: ", checkpoint_file)
      }

      cat("  Saved chunk", chunk, "of", R / chunk_size, "\n")
    }
  }

  cat("\nAll robustness checkpoints complete!\n")
}

# Save grid for reproducibility
write_csv(
  robustness_params,
  "results/extension3/robustness_simulation/parameter_grid.csv"
)

# Launch simulation
run_robustness(
  params = robustness_params,
  network_pools = new_pools_fast,
  R = 500,
  chunk_size = 100
)
