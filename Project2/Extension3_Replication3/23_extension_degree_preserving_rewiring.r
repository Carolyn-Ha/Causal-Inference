# ============================================================
# 23_extension_degree_preserving_rewiring.r
#
# Required input preparation for the final clustering analysis.
#
# Creates exactly the original a = 2 network pool used by script 24.
# Repeats the original rewiring draws (even though script 24 only
# uses the "original" member) to preserve the random-number
# sequence and therefore reproduce the originally analyzed pool.
#
# Output:
#   results/extension3/rewiring/matched_network_pools.rds
#
# Run this script before 24_extension_clustering_pairs.r.
# ============================================================

library(igraph)

set.seed(6165)

output_dir <- "results/extension3/rewiring"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_networks <- 150
n_nodes <- 50
attachment <- 2
m_edges <- 2
rewiring_multiplier <- 10

# Identical rewiring routine to the original script 23.
# The rewired graphs are retained in the output for reproducibility,
# but the final clustering analysis in script 24 uses only $original.
rewire_preserving_degrees <- function(g, multiplier = 10,
                                      max_attempts = 30) {
  for (attempt in seq_len(max_attempts)) {
    g_new <- rewire(
      g,
      with = keeping_degseq(
        niter = multiplier * ecount(g),
        loops = FALSE
      )
    )

    stopifnot(identical(degree(g), degree(g_new)))

    if (is_connected(g_new)) {
      return(g_new)
    }
  }

  stop("Could not generate a connected rewired network.")
}

original_networks <- vector("list", n_networks)
rewired_networks <- vector("list", n_networks)

# Keep this generation/rewiring interleaving unchanged:
# generating all original networks first would alter the RNG stream.
for (i in seq_len(n_networks)) {
  g_original <- sample_pa(
    n = n_nodes,
    power = attachment,
    m = m_edges,
    directed = FALSE,
    algorithm = "psumtree"
  )

  g_rewired <- rewire_preserving_degrees(
    g_original,
    multiplier = rewiring_multiplier
  )

  stopifnot(
    identical(degree(g_original), degree(g_rewired)),
    vcount(g_original) == vcount(g_rewired),
    ecount(g_original) == ecount(g_rewired)
  )

  original_networks[[i]] <- g_original
  rewired_networks[[i]] <- g_rewired

  if (i %% 25 == 0) {
    cat("Completed", i, "of", n_networks, "network pairs\n")
  }
}

saveRDS(
  list(original = original_networks, rewired = rewired_networks),
  file.path(output_dir, "matched_network_pools.rds")
)

cat("Saved original network pool for script 24.\n")
