library(fda)
library(MASS)
library(ggplot2)
library(splines)
library(parallel)

# import function_2.r

results_5 <- run_simulation_2(50, 5, iterations = 10)
plot(results_5)