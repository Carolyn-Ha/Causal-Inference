
library(igraph)


simulate_adoption <- function(g, seeds, alpha, beta, gamma, T = 3) {

  n <- vcount(g)

  # Adjacency matrix
  A <- as.matrix(as_adjacency_matrix(g))

  # Out-degree for calculating neighbor adoption
  d_out <- rowSums(A)
  A_norm <- A / pmax(d_out, 1)

  # In-degree for the degree-based seed effect
  d_in <- colSums(A)

  # Initialize adoption status
  adopted <- rep(0, n)
  adopted[seeds] <- 1

  # Total in-degree of selected seeds
  X <- sum(d_in[seeds])

  # Store adoption history
  adoption_history <- numeric(T + 1)
  adoption_history[1] <- mean(adopted)

  # Simulate adoption over T periods
  for (t in seq_len(T)) {

    # Fraction of neighbors who adopted previously
    Z <- as.vector(A_norm %*% adopted)

    # Random shocks
    epsilon <- rnorm(n)

    # Latent adoption propensity
    Y_star <- alpha + beta * Z + gamma * X + epsilon

    # Adoption is permanent
    adopted <- pmax(
      adopted,
      as.integer(Y_star > 0)
    )

    adoption_history[t + 1] <- mean(adopted)
  }

  return(list(
    final_adoption = mean(adopted),
    adoption_history = adoption_history,
    final_status = adopted
  ))
}
