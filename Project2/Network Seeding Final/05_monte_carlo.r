
library(tidyverse)
library(igraph)

source("01_seeding_functions.R")
source("02_adoption_model.R")
source("03_estimators.R")
source("04_simulate_experiment.R")




summarize_experiment <- function(df, rho = 0.5) {

  if (rho == 0) {

    estimates <- c(
      estimate_ht_offpolicy(df),
      estimate_hajek_offpolicy(df)
    )

    variances <- c(
      variance_ht_offpolicy(df),
      variance_hajek_offpolicy(df)
    )

    estimator_names <- c(
      "HT_offpolicy",
      "Hajek_offpolicy"
    )

  } else {

    estimates <- c(
      estimate_dm(df),
      estimate_ht(df),
      estimate_hajek(df)
    )

    variances <- c(
      variance_dm(df),
      variance_ht(df),
      variance_hajek(df)
    )

    estimator_names <- c(
      "DM",
      "HT",
      "Hajek"
    )
  }

  se <- sqrt(variances)
  z <- qnorm(0.95)

  tibble(
    estimator = estimator_names,
    estimate = estimates,
    variance = variances,
    SE = se,
    CI_lower = estimates - z * se,
    CI_upper = estimates + z * se,
    reject_H0 = ifelse(
      is.finite(se) & se > 0,
      abs(estimates / se) > z,
      NA
    )
  )
}






run_monte_carlo <- function(networks,
                            alpha,
                            beta,
                            gamma,
                            R = 100,
                            n_villages = length(networks),
                            rho = 0.5) {

  if (n_villages > length(networks)) {
    stop("n_villages cannot exceed the number of networks.")
  }

  results <- vector("list", R)

  for (r in seq_len(R)) {

    # Randomly select villages for this experiment
    selected <- sample(
      seq_along(networks),
      size = n_villages,
      replace = FALSE
    )

    experiment_networks <- networks[selected]

    # Run one experiment on the selected villages
    df_r <- simulate_experiment(
      networks = experiment_networks,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      rho = rho
    )

    # Calculate estimators and retain benchmark
    results[[r]] <- summarize_experiment(df_r, rho = rho) %>%
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




#set.seed(6165)

#mc_results <- run_monte_carlo(
  #networks = networks,
  #alpha = -2,
  #beta = 1.5,
  #gamma = 0.1,
  #R = 100
#)

# View the first few rows
#head(mc_results)

# Check the number of rows
#nrow(mc_results)



#power_results <- mc_results %>%
  #group_by(estimator) %>%
  #summarise(
    #n_reps = n(),
    #valid_tests = sum(!is.na(reject_H0)),
    #rejections = sum(reject_H0, na.rm = TRUE),
    #rejection_rate = if (
     # all(!is.na(reject_H0))
    #) {
     # mean(reject_H0)
    #} else {
     # NA_real_
    #},
    #mean_estimate = mean(estimate),
    #mean_SE = mean(SE),
    #.groups = "drop"
  #)

#print(power_results)



#ggplot(
 # power_results,
  #aes(
   # x = estimator,
    #y = rejection_rate,
    #fill = estimator
  #)
#) +
 # geom_col() +
  #scale_y_continuous(limits = c(0, 1)) +
  #labs(
   # title = "Rejection Rates Across 100 Simulated Experiments",
    #subtitle = "alpha = -2, beta = 1.5, gamma = 0.1",
    #x = "Estimator",
    #y = "Rejection rate"
  #) +
  #theme_minimal() +
  #theme(legend.position = "none")





#se_diagnostics <- mc_results %>%
 # group_by(estimator) %>%
  #summarise(
   # empirical_SD = sd(estimate),
    #mean_SE = mean(SE),
    #RMSE_SE_ratio = sqrt(mean(SE^2)) / sd(estimate),
    #.groups = "drop"
  #)

#print(se_diagnostics)


#df_off <- prepare_offpolicy(df)

#w <- df_off$p_A / df_off$p_B

#n_eff <- sum(w)^2 / sum(w^2)

#n_eff