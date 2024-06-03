library(fda)
library(MASS)
library(ggplot2)
library(splines)

# import function_3.r

nt_ny_list <- list(c(150, 10), c(150, 20), c(150, 5))
mse_results <- lapply(nt_ny_list, function(params) run_simulation(50, params[1], params[2], 10))