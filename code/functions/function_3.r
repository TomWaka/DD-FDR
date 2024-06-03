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
  
  for (i in 1:n) {
    xx[i, ] <- centered_gp_path(nn, 10, 10)
    xx[i, ] <- xx[i, ] - mean(xx[i, ])
    x[i, ] <- xx[i, xt] + rnorm(nt, 0, 1)
  }
  
  s <- 1:nn / nn
  t <- 1:nn / nn
  Beta <- matrix(nrow = nn, ncol = nn)
  for (t in 1:nn) {
    Beta[, t] <- GP_path(nn, 15, 10)
  }
  
  y <- matrix(nr = n, nc = ny)
  true <- matrix(nr = n, nc = nn)
  for (i in 1:n) {
    true[i, ] <- Beta %*% xx[i, ] / nn
    y[i, ] <- true[i, yt] + rnorm(ny, 0, 1)
  }
  
  list(x = x, y = y, true = true, xt = xt, yt = yt, Beta = Beta)
}

evaluate_model <- function(j, nbasis_x, nbasis_y, x, y, true, m_n1, xt, yt, nn) {
  Z <- matrix(nr = m_n1, nc = nbasis_x)
  V <- matrix(nr = m_n1, nc = nbasis_y)
  
  base_x <- ns(seq(0, 151, length = nn), df = nbasis_x)
  base_y <- ns(seq(0, 151, length = nn), df = nbasis_y)
  Phimat <- t(base_x) %*% (base_x)
  Psimat <- t(base_y) %*% (base_y)
  
  Bmat_x <- base_x[xt, ]
  Bmat_y <- base_y[yt, ]
  
  what <- ginv(t(Bmat_x) %*% Bmat_x) %*% t(Bmat_x) %*% t(x)
  Z <- t(what) %*% Phimat
  vhat <- ginv(t(Bmat_y) %*% Bmat_y) %*% t(Bmat_y) %*% t(y)
  V <- t(vhat)
  
  Ztr <- Z[1:m_n1, ]
  Vtr <- V[1:m_n1, ]
  Bhat <- ginv(Psimat %x% (t(Ztr) %*% Ztr)) %*% matrix(t(Ztr) %*% Vtr %*% Psimat, nc = 1)
  Bhat <- matrix(Bhat, nr = nbasis_x, nc = nbasis_y)
  
  yhat <- Z[-(1:m_n1), ] %*% Bhat %*% t(base_y)

  sqrt(mean(c(true[-(1:m_n1), ] - yhat)^2))
}

run_simulation <- function(m_n1, nt, ny, iterations = 10, nb = 24) {
  nn <- 150
  n <- 200
  mse <- numeric(nb)
  
  for (tt in 1:iterations) {
    print(tt)
    data <- generate_data(n, nn, nt, ny)
    x <- data$x
    y <- data$y
    true <- data$true
    xt <- data$xt
    yt <- data$yt
    
    results <- parallel::mclapply(1:nb, function(j) {
      nbasis_x <- 10
      nbasis_y <- 2 * j + 2
      evaluate_model(j, nbasis_x, nbasis_y, x, y, true, m_n1, xt, yt, nn)
    })
    
    for (j in 1:nb) {
      mse[j] <- mse[j] + results[[j]] / iterations
    }
  }
  mse
}

nt_ny_list <- list(c(150, 10), c(150, 20), c(150, 5))
mse_results <- lapply(nt_ny_list, function(params) run_simulation(50, params[1], params[2]))

# 結果のプロット
nbasis_values <- 2 + (1:24) * 2
data_plot <- data.frame(
  NB = rep(nbasis_values, 3),
  MSE = unlist(mse_results),
  district = factor(rep(c("M 10", "M 20", "M 5"), each = 24), levels = c("M 5", "M 10", "M 20"))
)

g <- ggplot(data_plot, aes(x = NB, y = MSE, color = district, group = district)) +
  geom_point(size = 1) +
  geom_line() +
  labs(x = "Number of bases of response", y = "MSE") + 
  coord_cartesian(ylim = c(0, 100)) +
  theme(
    legend.position = "bottom",
    legend.box.just = "right",
    legend.margin = margin(0, 0, 0, 0), legend.title = element_blank(),
    text = element_text(size = 15)
  )

# プロット表示と保存
plot(g)

pdf("FonF_1.pdf", width = 7, height = 4)
plot(g)
dev.off()
