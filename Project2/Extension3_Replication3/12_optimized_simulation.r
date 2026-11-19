# ============================================================
# 12_optimized_simulation.r
# Optimized reusable simulation functions (no automatic runs)
# ============================================================


# ==================================================
# 1. Precompute network information
# ==================================================

library(igraph)
library(tidyverse)

prepare_network <- function(g) {

  A <- as.matrix(as_adjacency_matrix(g))

  d_out <- rowSums(A)
  d_in <- colSums(A)

  # Normalized adjacency matrix for adoption
  A_norm <- A / pmax(d_out, 1)

  # One-hop selection probabilities
  valid_start <- d_out > 0

  A_valid <- A[valid_start, , drop = FALSE]
  d_valid <- d_out[valid_start]

  if (length(d_valid) == 0) {
    stop("No individuals have outgoing edges.")
  }

  q <- colSums(A_valid / d_valid) /
    nrow(A_valid)

  # Probability of selecting two distinct nodes
  pi_distinct <- 1 - sum(q^2)

  list(
    n = vcount(g),
    A_norm = A_norm,
    d_in = d_in,
    q = q,
    pi_distinct = pi_distinct,
    p_B = 1 / choose(vcount(g), 2)
  )
}



# ==================================================
# 2. Optimized adoption model
# ==================================================

simulate_adoption_fast <- function(
    net, seeds, alpha, beta, gamma, T = 3
) {

  adopted <- numeric(net$n)
  adopted[seeds] <- 1

  X <- sum(net$d_in[seeds])

  for (t in seq_len(T)) {

    Z <- as.vector(net$A_norm %*% adopted)

    epsilon <- rnorm(net$n)

    Y_star <- alpha + beta * Z + gamma * X + epsilon

    adopted <- pmax(
      adopted,
      as.integer(Y_star > 0)
    )
  }

  mean(adopted)
}



# ==================================================
# 3. Optimized experiment function
# ==================================================

simulate_experiment_fast <- function(
    networks,
    alpha,
    beta,
    gamma,
    rho = 0.5,
    T = 3
) {

  N <- length(networks)

  strategy <- character(N)
  outcome <- numeric(N)
  y_A <- numeric(N)
  y_B <- numeric(N)
  p_A <- numeric(N)
  p_B <- numeric(N)
  p_design <- numeric(N)

  for (i in seq_len(N)) {

    net <- networks[[i]]
    n <- net$n

    # Assign policy
    strategy[i] <- ifelse(
      rbinom(1, 1, rho) == 1,
      "A", "B"
    )

    # Draw one-hop seed set
    seeds_A <- sample(
      seq_len(n),
      size = 2,
      replace = TRUE,
      prob = net$q
    )

    while (length(unique(seeds_A)) < 2) {
      seeds_A <- sample(
        seq_len(n),
        size = 2,
        replace = TRUE,
        prob = net$q
      )
    }

    # Draw random seed set
    seeds_B <- sample(
      seq_len(n),
      size = 2,
      replace = FALSE
    )

    # Actual observed seed set
    seeds_observed <- if (strategy[i] == "A") {
      seeds_A
    } else {
      seeds_B
    }

    # Seed-set probabilities
    p_A[i] <- 2 * prod(net$q[seeds_observed]) /
      net$pi_distinct

    p_B[i] <- net$p_B

    p_design[i] <- rho * p_A[i] +
      (1 - rho) * p_B[i]

    # Adoption under both strategies
    y_A[i] <- simulate_adoption_fast(
      net, seeds_A, alpha, beta, gamma, T
    )

    y_B[i] <- simulate_adoption_fast(
      net, seeds_B, alpha, beta, gamma, T
    )

    # Observed outcome
    outcome[i] <- if (strategy[i] == "A") {
      y_A[i]
    } else {
      y_B[i]
    }
  }

  tibble(
    village_id = seq_len(N),
    strategy = strategy,
    outcome = outcome,
    y_A = y_A,
    y_B = y_B,
    p_A = p_A,
    p_B = p_B,
    p_design = p_design
  )
}



# ==================================================
# 4. Optimized Monte Carlo runner
# ==================================================

run_monte_carlo_fast <- function(
    networks,
    alpha,
    beta,
    gamma,
    R = 100,
    n_villages = length(networks),
    rho = 0.5
) {

  if (n_villages > length(networks)) {
    stop("Not enough networks.")
  }

  results <- vector("list", R)

  for (r in seq_len(R)) {

    # Randomly select villages
    selected <- sample(
      seq_along(networks),
      size = n_villages,
      replace = FALSE
    )

    experiment_networks <- networks[selected]

    df_r <- simulate_experiment_fast(
      networks = experiment_networks,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      rho = rho
    )

    results[[r]] <- summarize_experiment(
      df_r,
      rho = rho
    ) %>%
      mutate(
        rep = r,
        alpha = alpha,
        beta = beta,
        gamma = gamma,
        n_villages = n_villages,
        rho = rho,
        tau_benchmark = mean(df_r$y_A - df_r$y_B)
      )
  }

  bind_rows(results)
}
