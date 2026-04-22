# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  Dinosaur GLM Phylodynamics
#        V\ Y /V    Get priors 95% interval
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# Load libraries and functions -------------------------------------------------
library(tidyverse)

get_prior_table <- function(param, distr, sd, meanl, sdl, upper, lower, mean, value, 
                            n_round, n_samples = 25000) {
  if (distr == "LogNormal") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", meanl, ", ", sdl, ")"),
                median = round(qlnorm(0.5, meanl, sdl), n_round),
                density95 = paste0("(", round(qlnorm(0.025, meanl, sdl), n_round)," - ", 
                                   round(qlnorm(0.975, meanl, sdl), n_round), ")"),
                samples =  list(rlnorm(n_samples, meanlog = meanl, sdlog = sdl))
    )
  } else if (distr == "Normal") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", mean, ", ", sd, ")"),
                median = round(qnorm(0.5, mean, sd), n_round),
                density95 = paste0("(", round(qnorm(0.025, mean, sd), n_round)," - ", 
                                   round(qnorm(0.975, mean, sd), n_round), ")"),
                samples = list(rnorm(n_samples, mean, sd))
    )
  } else if (distr == "Uniform") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", lower, ", ", upper, ")"),
                median = round(qunif(0.5, lower, upper), n_round),
                density95 = paste0("(", round(qunif(0.025, lower, upper), n_round)," - ", 
                                   round(qunif(0.975, lower, upper), n_round), ")"),
                samples = list(runif(n_samples, min = lower, max = upper))
    )
  } else if (distr == "OneOnX") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", lower, ", ", upper, ")"),
                #median = round(qunif(0.5, lower, upper), n_round),
                #density95 = paste0("(", round(qunif(0.025, lower, upper), n_round)," - ", 
                #                   round(qunif(0.975, lower, upper), n_round), ")")
    )
  } else if (distr == "Exponential") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, "(", mean, ")"),
                median = round(qexp(0.5, mean), n_round),
                density95 = paste0("(", round(qexp(0.025, mean), n_round)," - ", 
                                   round(qexp(0.975, mean), n_round), ")"),
                samples = list(rexp(n_samples, rate = mean))
    )
  } else if (distr == "Fixed") {
    t <- tibble(parameter = param,
                distribution = paste0(distr, " to ", value))
  } else {
    t <- tibble(parameter = param,
                distribution = distr)
  }
  return(t)
}


# Priors Dinosaur analyses -----------------------------------------------------

priors <- bind_rows(
  get_prior_table("Origin", "Uniform", lower = 0, upper = 200, n_round = 3), 
  get_prior_table("Death rate", "Exponential", mean = 1.0, n_round = 3), 
  get_prior_table("Birth rate", "Exponential", mean = 1.0, n_round = 3), 
  get_prior_table("Sampling rate", "Exponential", mean = 0.2, n_round = 3),
  get_prior_table("GLM Baseline", "Exponential", mean = 0.2, n_round = 3), 
  get_prior_table("GLM Coefficients", "Normal", mean = 0.0, sd =  1.0, n_round = 3), 
  get_prior_table("Indicator Sum", "Poisson", mean = 0.693, n_round = 3))
# 
# print(xtable::xtable(priors %>% select(-samples), type = "latex"), file = "results/report/priors.tex")

