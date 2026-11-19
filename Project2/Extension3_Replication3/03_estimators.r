
# ============================================================
# 03_estimators.r
# DM, HT, and Hajek estimators and variance estimators
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Difference-in-means estimator
# ------------------------------------------------------------

estimate_dm <- function(df) {

  mean_A <- mean(df$outcome[df$strategy == "A"])
  mean_B <- mean(df$outcome[df$strategy == "B"])

  return(mean_A - mean_B)
}


# ------------------------------------------------------------
# 2. Horvitz-Thompson estimator
# ------------------------------------------------------------

estimate_ht <- function(df) {

  w_A <- df$p_A / df$p_design
  w_B <- df$p_B / df$p_design

  mean((w_A - w_B) * df$outcome)
}


# ------------------------------------------------------------
# 3. Hajek estimator
# ------------------------------------------------------------

estimate_hajek <- function(df) {

  w_A <- df$p_A / df$p_design
  w_B <- df$p_B / df$p_design

  mean_A <- sum(w_A * df$outcome) / sum(w_A)
  mean_B <- sum(w_B * df$outcome) / sum(w_B)

  return(mean_A - mean_B)
}


# ------------------------------------------------------------
# 4. Horvitz-Thompson variance estimator
# ------------------------------------------------------------

variance_ht <- function(df) {

  w_A <- df$p_A / df$p_design
  w_B <- df$p_B / df$p_design

  tau_hat <- estimate_ht(df)

  influence <- (w_A - w_B) * df$outcome

  V_hat <- mean((influence - tau_hat)^2)

  return(V_hat / nrow(df))
}


# ------------------------------------------------------------
# 5. Difference-in-means variance estimator
# ------------------------------------------------------------

variance_dm <- function(df) {

  y_A <- df$outcome[df$strategy == "A"]
  y_B <- df$outcome[df$strategy == "B"]

  n_A <- length(y_A)
  n_B <- length(y_B)

  if (n_A < 2 || n_B < 2) {
    return(NA_real_)
  }

  var(y_A) / n_A + var(y_B) / n_B
}


# ------------------------------------------------------------
# 6. Hajek variance estimator
# ------------------------------------------------------------

variance_hajek <- function(df) {

  n <- nrow(df)

  w_A <- df$p_A / df$p_design
  w_B <- df$p_B / df$p_design

  mu_A <- sum(w_A * df$outcome) / sum(w_A)
  mu_B <- sum(w_B * df$outcome) / sum(w_B)

  influence <- w_A * (df$outcome - mu_A) -
               w_B * (df$outcome - mu_B)

  V_hat <- mean(influence^2)

  return(V_hat / n)
}


# ------------------------------------------------------------
# 7. Prepare off-policy dataset
# ------------------------------------------------------------

prepare_offpolicy <- function(df) {

  df_off <- df %>%
    filter(strategy == "B") %>%
    mutate(p_design = p_B)

  return(df_off)
}


# ------------------------------------------------------------
# 8. Off-policy Horvitz-Thompson estimator
# ------------------------------------------------------------

estimate_ht_offpolicy <- function(df) {

  df_off <- prepare_offpolicy(df)

  estimate_ht(df_off)
}


# ------------------------------------------------------------
# 9. Off-policy Hajek estimator
# ------------------------------------------------------------

estimate_hajek_offpolicy <- function(df) {

  df_off <- prepare_offpolicy(df)

  estimate_hajek(df_off)
}


# ------------------------------------------------------------
# 10. Off-policy Horvitz-Thompson variance
# ------------------------------------------------------------

variance_ht_offpolicy <- function(df) {

  df_off <- prepare_offpolicy(df)

  variance_ht(df_off)
}


# ------------------------------------------------------------
# 11. Off-policy Hajek variance
# ------------------------------------------------------------

variance_hajek_offpolicy <- function(df) {

  df_off <- prepare_offpolicy(df)

  variance_hajek(df_off)
}

