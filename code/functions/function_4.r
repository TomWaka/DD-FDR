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

generate_data <- function(n, nn, nt, ny) {
  xt <- sort(sample(1:nn, nt))
  yt <- sort(sample(1:nn, ny))
  xx <- matrix(0, nrow = n, ncol = nn)
  x <- matrix(0, nrow = n, ncol = nt)
  y <- matrix(0, nrow = n, ncol = ny)
  
  for (i in 1:n) {
    xx[i, ] <- centered_gp_path(nn, 10, 10)
    x[i, ] <- xx[i, xt] + rnorm(nt, 0, 1)
  }
  
  Beta <- matrix(nrow = nn, ncol = nn)
  for (t in 1:nn) {
    Beta[, t] <- GP_path(nn, 15, 10)
  }
  
  true <- matrix(nrow = n, ncol = nn)
  for (i in 1:n) {
    true[i, ] <- Beta %*% xx[i, ] / nn
    y[i, ] <- true[i, yt] + rnorm(ny, 0, 1)
  }
  
  list(x = x, y = y, true = true, xt = xt, yt = yt, Beta = Beta)
}

evaluate_model <- function(nbasis_x, nbasis_y, x, y, true, xt, yt, nn, m_n1) {
  Z <- matrix(nrow = m_n1, ncol = nbasis_x)
  V <- matrix(nrow = m_n1, ncol = nbasis_y)
  
  base_x <- ns(seq(0, 151, length = nn), df = nbasis_x)
  base_y <- ns(seq(0, 151, length = nn), df = nbasis_y)
  Phimat <- crossprod(base_x)
  Psimat <- crossprod(base_y)

  Bmat_x <- base_x[xt, ]
  Bmat_y <- base_y[yt, ]
  
  what <- ginv(crossprod(Bmat_x)) %*% t(Bmat_x) %*% t(x)
  Z <- t(what) %*% Phimat
  
  vhat <- ginv(crossprod(Bmat_y)) %*% t(Bmat_y) %*% t(y)
  V <- t(vhat)
  
  Ztr <- Z[1:m_n1, ]
  Vtr <- V[1:m_n1, ]

  RR <- Psimat %x% (t(Ztr) %*% Ztr)
  Bhat <- ginv( RR + 5e-6 * diag(nrow(RR))) %*% matrix(t(Ztr) %*% Vtr %*% Psimat, nc = 1)
  Bhat <- matrix(Bhat, nrow = nbasis_x, ncol = nbasis_y)
  
  yhat <- Z[-(1:m_n1), ] %*% Bhat %*% t(base_y)

  sqrt(mean((true[-(1:m_n1), ] - yhat)^2))
}

run_simulation <- function(m_n1, nt, ny, nb = 24, iterations = 10) {
  nn <- 150
  n <- 200
  mse <- numeric(nb)
  num <- 0
  
  for (tt in 1:iterations) {
    data <- generate_data(n, nn, nt, ny)
    x <- data$x
    y <- data$y
    true <- data$true
    xt <- data$xt
    yt <- data$yt
    
    results <- parallel::mclapply(1:nb, function(j) {
      nbasis_x <- 2 * j + 2
      nbasis_y <- 10
      evaluate_model(nbasis_x, nbasis_y, x, y, true, xt, yt, nn, m_n1)
    })

    for (j in 1:nb) {
      mse[j] <- mse[j] + results[[j]] / iterations
    }
  }
  mse
}