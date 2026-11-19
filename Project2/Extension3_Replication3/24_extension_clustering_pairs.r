
# ============================================================
# 24_extension_clustering_pairs.r
#
# Goal:
# Create matched low-clustering and high-clustering
# networks while preserving every individual node's degree.
#
# Each pair comes from the SAME original network.
#
# Outputs:
#   - Matched network pools
#   - Network metrics
#   - Paired differences
#   - Diagnostic plots
#
# No adoption simulations are performed in this script.
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

library(igraph)
library(tidyverse)

set.seed(6165)

input_file <- paste0(
  "results/extension3/rewiring/",
  "matched_network_pools.rds"
)

output_dir <- "results/extension3/clustering_pairs"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Number of independently rewired candidates per graph
n_candidates <- 60

# Number of attempted swaps per candidate, as a
# multiple of the original graph's edge count
rewiring_multiplier <- 10


# ------------------------------------------------------------
# 2. Load our existing original networks
# ------------------------------------------------------------

saved_networks <- readRDS(input_file)

# Reuse the original graphs, not their earlier rewired versions
original_networks <- saved_networks$original

n_networks <- length(original_networks)

cat("Loaded", n_networks, "original networks\n")

stopifnot(n_networks == 150)


# ------------------------------------------------------------
# 3. Function to calculate network metrics
# ------------------------------------------------------------

network_metrics <- function(g) {

  tibble(
    degree_sd = sd(degree(g)),

    clustering = transitivity(
      g,
      type = "global"
    ),

    avg_path_length = mean_distance(
      g,
      directed = FALSE,
      unconnected = TRUE
    ),

    assortativity = assortativity_degree(
      g,
      directed = FALSE
    ),

    connected = is_connected(g)
  )
}


# ------------------------------------------------------------
# 4. Generate high- and low-clustering versions
# ------------------------------------------------------------

create_clustering_pair <- function(
    g,
    n_candidates = 60,
    rewiring_multiplier = 10
) {

  stopifnot(
    !is_directed(g),
    is_simple(g),
    is_connected(g)
  )

  original_degrees <- degree(g)

  # Include the original as an additional candidate
  candidates <- list(g)

  # Generate independent rewired candidates
  for (j in seq_len(n_candidates)) {

    candidate <- rewire(
      g,
      with = keeping_degseq(
        niter = rewiring_multiplier * ecount(g),
        loops = FALSE
      )
    )

    # Check degree preservation
    stopifnot(
      identical(
        degree(candidate),
        original_degrees
      )
    )

    # Retain only connected, simple graphs
    if (
      is_connected(candidate) &&
      is_simple(candidate)
    ) {

      candidates[[length(candidates) + 1]] <- candidate
    }
  }

  # Measure clustering of every accepted candidate
  clustering_values <- vapply(
    candidates,
    function(net) {
      transitivity(net, type = "global")
    },
    numeric(1)
  )

  # Remove undefined clustering values, if any
  valid <- which(is.finite(clustering_values))

  if (length(valid) < 2) {
    stop("Fewer than two valid clustering candidates.")
  }

  # Select extremes among available candidates
  low_index <- valid[
    which.min(clustering_values[valid])
  ]

  high_index <- valid[
    which.max(clustering_values[valid])
  ]

  g_low <- candidates[[low_index]]
  g_high <- candidates[[high_index]]

  # Verify degree sequences match exactly
  stopifnot(
    identical(degree(g_low), degree(g_high)),
    identical(degree(g_low), original_degrees)
  )

  list(
    low = g_low,
    high = g_high,
    n_valid_candidates = length(valid),
    low_clustering = clustering_values[low_index],
    high_clustering = clustering_values[high_index],
    clustering_gap =
      clustering_values[high_index] -
      clustering_values[low_index]
  )
}


# ------------------------------------------------------------
# 5. Create 150 matched network pairs
# ------------------------------------------------------------

low_networks <- vector("list", n_networks)
high_networks <- vector("list", n_networks)

pair_info <- vector("list", n_networks)

for (i in seq_len(n_networks)) {

  pair <- create_clustering_pair(
    original_networks[[i]],
    n_candidates = n_candidates,
    rewiring_multiplier = rewiring_multiplier
  )

  low_networks[[i]] <- pair$low
  high_networks[[i]] <- pair$high

  pair_info[[i]] <- tibble(
    network_id = i,
    n_valid_candidates = pair$n_valid_candidates,
    clustering_low = pair$low_clustering,
    clustering_high = pair$high_clustering,
    clustering_gap = pair$clustering_gap
  )

  if (i %% 10 == 0) {
    cat(
      "Completed", i, "of", n_networks,
      "matched pairs\n"
    )
  }
}

pair_info <- bind_rows(pair_info)

# Check whether all pairs have positive clustering gaps
cat(
  "\nPairs with positive clustering difference:",
  sum(pair_info$clustering_gap > 0),
  "out of", n_networks, "\n"
)

if (any(pair_info$clustering_gap <= 0)) {
  warning(
    "Some pairs have identical clustering. ",
    "Consider increasing n_candidates."
  )
}


# ------------------------------------------------------------
# 6. Calculate metrics for both network groups
# ------------------------------------------------------------

low_metrics <- map2_dfr(
  low_networks,
  seq_len(n_networks),
  function(g, i) {

    network_metrics(g) %>%
      mutate(
        network_id = i,
        condition = "Low clustering"
      )
  }
)

high_metrics <- map2_dfr(
  high_networks,
  seq_len(n_networks),
  function(g, i) {

    network_metrics(g) %>%
      mutate(
        network_id = i,
        condition = "High clustering"
      )
  }
)

all_metrics <- bind_rows(
  low_metrics,
  high_metrics
)


# ------------------------------------------------------------
# 7. Calculate differences within matched pairs
# ------------------------------------------------------------

paired_metrics <- low_metrics %>%
  select(
    network_id,
    degree_sd,
    clustering,
    avg_path_length,
    assortativity
  ) %>%
  inner_join(
    high_metrics %>%
      select(
        network_id,
        degree_sd,
        clustering,
        avg_path_length,
        assortativity
      ),
    by = "network_id",
    suffix = c("_low", "_high")
  ) %>%
  mutate(
    delta_degree_sd =
      degree_sd_high - degree_sd_low,

    delta_clustering =
      clustering_high - clustering_low,

    delta_path_length =
      avg_path_length_high - avg_path_length_low,

    delta_assortativity =
      assortativity_high - assortativity_low
  )

# Confirm degree heterogeneity is unchanged
stopifnot(
  all(paired_metrics$delta_degree_sd == 0)
)

# Confirm clustering never decreases
stopifnot(
  all(paired_metrics$delta_clustering >= 0)
)


# ------------------------------------------------------------
# 8. Summarize network characteristics
# ------------------------------------------------------------

metrics_summary <- all_metrics %>%
  group_by(condition) %>%
  summarise(
    mean_degree_sd = mean(degree_sd),

    mean_clustering = mean(
      clustering,
      na.rm = TRUE
    ),

    mean_path_length = mean(
      avg_path_length,
      na.rm = TRUE
    ),

    mean_assortativity = mean(
      assortativity,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

cat("\nAverage network characteristics:\n")
print(metrics_summary, width = Inf)


# ------------------------------------------------------------
# 9. Summarize paired differences
# ------------------------------------------------------------

change_summary <- paired_metrics %>%
  select(
    network_id,
    starts_with("delta_")
  ) %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "metric",
    values_to = "change"
  ) %>%
  group_by(metric) %>%
  summarise(
    mean_change = mean(change, na.rm = TRUE),

    mean_absolute_change = mean(
      abs(change),
      na.rm = TRUE
    ),

    median_change = median(change, na.rm = TRUE),

    sd_change = sd(change, na.rm = TRUE),

    proportion_increased = mean(
      change > 0,
      na.rm = TRUE
    ),

    min_change = min(change, na.rm = TRUE),

    max_change = max(change, na.rm = TRUE),

    .groups = "drop"
  )

cat("\nSummary of matched changes:\n")
print(change_summary, width = Inf)


# ------------------------------------------------------------
# 10. Plot original versus high/low network metrics
# ------------------------------------------------------------

metrics_long <- all_metrics %>%
  select(
    network_id,
    condition,
    degree_sd,
    clustering,
    avg_path_length,
    assortativity
  ) %>%
  pivot_longer(
    cols = c(
      degree_sd,
      clustering,
      avg_path_length,
      assortativity
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      degree_sd = "Degree SD",
      clustering = "Global clustering",
      avg_path_length = "Average path length",
      assortativity = "Degree assortativity"
    )
  )

p_metrics <- ggplot(
  metrics_long,
  aes(
    x = condition,
    y = value,
    fill = condition
  )
) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Degree-Preserving Clustering Comparison",
    subtitle = paste(
      n_networks,
      "matched high- and low-clustering pairs"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p_metrics)

ggsave(
  file.path(output_dir, "clustering_groups_metrics.png"),
  p_metrics,
  width = 10,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------
# 11. Plot paired changes in structural characteristics
# ------------------------------------------------------------

changes_long <- paired_metrics %>%
  select(
    network_id,
    delta_clustering,
    delta_path_length,
    delta_assortativity
  ) %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "metric",
    values_to = "change"
  )

p_changes <- ggplot(
  changes_long,
  aes(x = change)
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
  facet_wrap(
    ~ metric,
    scales = "free_x",
    ncol = 1
  ) +
  labs(
    title = "Structural Differences Within Matched Pairs",
    subtitle = "High clustering minus low clustering",
    x = "Paired difference",
    y = "Number of network pairs"
  ) +
  theme_minimal(base_size = 12)

print(p_changes)

ggsave(
  file.path(output_dir, "clustering_paired_changes.png"),
  p_changes,
  width = 9,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 12. Save everything for our adoption simulations
# ------------------------------------------------------------

saveRDS(
  list(
    original = original_networks,
    low = low_networks,
    high = high_networks,
    pair_info = pair_info
  ),
  file.path(output_dir, "clustering_network_pools.rds")
)

write_csv(
  all_metrics,
  file.path(output_dir, "clustering_network_metrics.csv")
)

write_csv(
  paired_metrics,
  file.path(output_dir, "clustering_paired_metrics.csv")
)

write_csv(
  change_summary,
  file.path(output_dir, "clustering_change_summary.csv")
)

write_csv(
  pair_info,
  file.path(output_dir, "clustering_pair_info.csv")
)

cat("\nClustering-pair generation complete.\n")
cat("Outputs saved to:", output_dir, "\n")
