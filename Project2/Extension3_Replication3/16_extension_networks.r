
# ==================================================
# 16_extension_networks.r
# Extension 3: Generate and prepare synthetic network pools
# ==================================================

library(tidyverse)
library(igraph)

# Functions for preparing each generated graph for the fast simulator.
source("12_optimized_simulation.r")

set.seed(6165)

# Network settings
n_nodes <- 50
m_edges <- 2

# Levels of preferential attachment
attachment_levels <- c(0, 0.5, 1, 1.5, 2)

# Generate one network at each level
example_networks <- lapply(
  attachment_levels,
  function(a) {
    sample_pa(
      n = n_nodes,
      power = a,
      m = m_edges,
      directed = FALSE,
      algorithm = "psumtree"
    )
  }
)

names(example_networks) <- paste0(
  "attachment_", attachment_levels
)

# Measure the resulting network structure
network_summary <- tibble(
  attachment = attachment_levels,
  nodes = sapply(example_networks, vcount),
  edges = sapply(example_networks, ecount),
  mean_degree = sapply(
    example_networks,
    function(g) mean(degree(g))
  ),
  max_degree = sapply(
    example_networks,
    function(g) max(degree(g))
  ),
  degree_sd = sapply(
    example_networks,
    function(g) sd(degree(g))
  ),
  clustering = sapply(
    example_networks,
    function(g) transitivity(g, type = "global")
  )
)

print(network_summary, width = Inf)




# ==================================================
# Step 2. Generate synthetic village pools
# ==================================================

set.seed(6166)

n_villages_pool <- 150

# Generate 150 networks for each attachment level
network_pools <- lapply(
  attachment_levels,
  function(a) {

    lapply(
      seq_len(n_villages_pool),
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

names(network_pools) <- paste0(
  "attachment_", attachment_levels
)

# Summarize average structural characteristics
pool_summary <- purrr::map2_dfr(
  network_pools,
  attachment_levels,
  function(pool, a) {

    tibble(
      attachment = a,
      mean_degree_sd = mean(
        sapply(pool, function(g) sd(degree(g)))
      ),
      mean_max_degree = mean(
        sapply(pool, function(g) max(degree(g)))
      ),
      mean_clustering = mean(
        sapply(pool, function(g) {
          transitivity(g, type = "global")
        }),
        na.rm = TRUE
      )
    )
  }
)

print(pool_summary, width = Inf)





# ==================================================
# Step 3. Validate generated network structures
# ==================================================

network_checks <- purrr::map2_dfr(
  network_pools,
  attachment_levels,
  function(pool, a) {
    tibble(
      attachment = a,
      n_networks = length(pool),
      non_simple = sum(
        !sapply(pool, igraph::is_simple)
      ),
      disconnected = sum(
        !sapply(pool, igraph::is_connected)
      )
    )
  }
)

print(network_checks)



# ==================================================
# Step 4. Prepare and save synthetic network pools
# ==================================================

# Use the optimized prepare_network() function
# from 12_optimized_simulation.r

network_pools_fast <- lapply(
  network_pools,
  function(pool) {
    lapply(pool, prepare_network)
  }
)

# Save prepared networks
dir.create(
  "results/extension3",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  network_pools_fast,
  "results/extension3/network_pools_fast.rds"
)

# Verify one-hop sampling probabilities
probability_checks <- purrr::map_dfr(
  network_pools_fast,
  function(pool) {
    tibble(
      min_prob_sum = min(sapply(pool, function(net) sum(net$q))),
      max_prob_sum = max(sapply(pool, function(net) sum(net$q))),
      min_distinct_prob = min(sapply(pool, function(net) net$pi_distinct))
    )
  },
  .id = "network_type"
)

print(probability_checks)

