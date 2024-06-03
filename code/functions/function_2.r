library(fda)
library(MASS)
library(ggplot2)
library(splines)
library(parallel)

GP_path <- function(K, th, h) {
  grid <- seq(1, K)
  dist_matrix <- outer(grid, grid, FUN = function(x1, x2) (x1 - x2)^2)
  gram <- th^2 * exp(-dist_matrix / (2 * h^2))
  mvrnorm(1, mu = rep(0, K), Sigma = gram)
}

centered_gp_path <- function(nn, lscale, var) {
  path <- GP_path(nn, lscale, var)
  path - mean(path)
}

generate_data <- function(n, nn, nt) {
  xt <- sort(sample(1:nn, nt))
  xx <- matrix(0, nrow = n, ncol = nn)
  x <- matrix(0, nrow = n, ncol = nt)
  
  for (i in 1:n) {
    xx[i, ] <- centered_gp_path(nn, 10, 10)
    x[i, ] <- xx[i, xt] + rnorm(nt, 0, 1)
  }
  
  bb <- GP_path(nn, 15, 10)
  true <- apply(xx, 1, function(row) sum(bb * row) / nn)
  y <- true + rnorm(n, 0, 1)
  
  list(x = x, y = y, true = true, xt = xt)
}

evaluate_model <- function(nbasis, x, y, true, m_n1, xt, nn, nt) {
  Z <- matrix(nrow = nrow(x), ncol = nbasis)
  BBmat <- ns(seq(0, 151, length = nn), df = nbasis)
  Bmat <- BBmat[xt, ]
  Phimat <- crossprod(BBmat)
  giB <- ginv(crossprod(Bmat))
  
  for (i in 1:nrow(x)) {
    what <- giB %*% t(Bmat) %*% x[i, ]
    Z[i, ] <- t(what) %*% Phimat
  }
  
  bhat <- ginv(crossprod(Z[1:m_n1, ])) %*% crossprod(Z[1:m_n1, ], y[1:m_n1])
  yhat <- Z[-(1:m_n1), ] %*% bhat
  sqrt(mean((true[-(1:m_n1)] - yhat)^2))
}

run_simulation_2 <- function(m_n1, nt, iterations = 100, m_M0 = 24) {
  nn <- 150
  n <- 200
  mse <- numeric(m_M0)
  
  for (iteration in 1:iterations) {
    data <- generate_data(n, nn, nt)
    x <- data$x
    y <- data$y
    true <- data$true
    xt <- data$xt
    
    results <- parallel::mclapply(1:m_M0, function(j) {
      nbasis <- 2 * j + 2
      result <- evaluate_model(nbasis, x, y, true, m_n1, xt, nn, nt)
      result
    })
    
    for (j in 1:m_M0) {
      mse[j] <- mse[j] + results[[j]] / iterations
    }
  }
  mse
}
