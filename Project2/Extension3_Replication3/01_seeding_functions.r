# ============================================================
# 01_seeding_functions.r
# Calculate one-hop seed selection probabilities
# ============================================================

one_hop_probs <- function(g) {

  A <- as.matrix(igraph::as_adjacency_matrix(g))

  d_out <- rowSums(A)

  valid_start <- d_out > 0

  A_valid <- A[valid_start, , drop = FALSE]
  d_valid <- d_out[valid_start]

  if (length(d_valid) == 0) {
    stop("No individuals have outgoing edges.")
  }

  q <- colSums(A_valid / d_valid) / nrow(A_valid)

  return(q)
}