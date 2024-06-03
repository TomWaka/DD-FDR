library(fda)
library(MASS)
library(ggplot2)
library(splines)

# import function_4.r

nt_ny_list <- c(5, 10, 20)
mse_results <- lapply(nt_ny_list, function(m_n1) run_simulation(m_n1, 150, 75, nb = 24, iterations = 10))
