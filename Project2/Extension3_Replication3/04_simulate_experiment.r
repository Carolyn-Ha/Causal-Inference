
# ============================================================
# 04_simulate_experiment.r
# Simulate experiments comparing one-hop and random seeding
# ============================================================

library(igraph)
library(tidyverse)

# ------------------------------------------------------------
# 1. Simulate one experiment across a collection of networks
# ------------------------------------------------------------

simulate_experiment <- function(networks,
                                alpha,
                                beta,
                                gamma,
                                k = 2,
                                rho = 0.5,
                                T = 3) {

  if (k != 2) {
    stop("This version only supports k = 2.")
  }

  results <- vector("list", length(networks))

  for (i in seq_along(networks)) {

    g <- networks[[i]]
    n <- vcount(g)

    # 1. Assign village to strategy A or B
    strategy <- ifelse(
      rbinom(1, 1, rho) == 1,
      "A", "B"
    )

    # 2. Calculate single-node one-hop probabilities
    q <- one_hop_probs(g)

    # 3. Draw a seed set under strategy A
    seeds_A <- sample(
      seq_len(n),
      size = k,
      replace = TRUE,
      prob = q
    )

    while (length(unique(seeds_A)) < k) {
      seeds_A <- sample(
        seq_len(n),
        size = k,
        replace = TRUE,
        prob = q
      )
    }

    # 4. Draw a seed set under strategy B
    seeds_B <- sample(
      seq_len(n),
      size = k,
      replace = FALSE
    )

    # 5. Select which seed set is observed
    seeds_observed <- if (strategy == "A") {
      seeds_A
    } else {
      seeds_B
    }

    # 6. Calculate probabilities of the observed seed set
    p_A <- 2 * prod(q[seeds_observed]) /
      (1 - sum(q^2))

    p_B <- 1 / choose(n, k)

    p_design <- rho * p_A + (1 - rho) * p_B

    # 7. Simulate adoption under strategy A
    adoption_A <- simulate_adoption(
      g = g,
      seeds = seeds_A,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      T = T
    )

    # 8. Simulate adoption under strategy B
    adoption_B <- simulate_adoption(
      g = g,
      seeds = seeds_B,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      T = T
    )

    y_A <- adoption_A$final_adoption
    y_B <- adoption_B$final_adoption

    # 9. Record only the assigned outcome
    outcome <- if (strategy == "A") y_A else y_B

    # 10. Save village-level information
    results[[i]] <- tibble(
      village_id = i,
      strategy = strategy,
      outcome = outcome,
      y_A = y_A,
      y_B = y_B,
      p_A = p_A,
      p_B = p_B,
      p_design = p_design
    )
  }

  bind_rows(results)
}


