packages <- c("refund", "fda.usc", "ggplot2", "splines")

installed_packages <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

evaluate_model <- function(x, y, m_n1, n, m_M0, nt, nn) {
  MSE <- numeric(m_M0)
  
  for (j in 1:m_M0) {
    nbasis <- j + 4
    Z <- matrix(nr = n, nc = nbasis)
    BBmat <- ns(seq(0, 401, length = nn), df = nbasis)
    Bmat <- BBmat[1:nt, ]
    Phimat <- crossprod(BBmat)
    
    giB <- ginv(crossprod(Bmat))
    for (i in 1:n) {
      what <- giB %*% t(Bmat) %*% x[i, ]
      Z[i, ] <- t(what) %*% Phimat
    }
    
    bhat <- ginv(crossprod(Z[1:m_n1, ])) %*% t(Z[1:m_n1, ]) %*% y[1:m_n1]
    yhat <- Z[-(1:m_n1), ] %*% bhat
    MSE[j] <- sqrt(mean((y[-(1:m_n1)] - yhat)^2))
  }
  
  return(MSE)
}

evaluate_and_plot <- function(x, y, m_n1, n, m_M0, nt, nn, file_name) {
  MSE <- evaluate_model(x, y, m_n1, n, m_M0, nt, nn)
  
  plot_data <- data.frame(Nbasis = 1:m_M0 + 4, MSE = MSE)
  
  p <- ggplot(plot_data, aes(x = Nbasis, y = MSE)) +
    geom_point() +
    geom_line() +
    labs(x = "Number of Bases", y = "MSE") +
    theme(
      legend.position = "bottom",
      legend.box.just = "right",
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank(),
      text = element_text(size = 15)
    )
  
  print(p)
  
  pdf(file_name, width = 7, height = 4)
  print(p)
  dev.off()
}

# gasoline
data(gasoline)
x_gasoline <- gasoline$NIR
y_gasoline <- gasoline$octane
evaluate_and_plot(x_gasoline, y_gasoline, m_n1 = 10, n = 60, m_M0 = 60, nt = 401, nn = 401, file_name = "gasoline.pdf")

# DTI
data(DTI)
y_DTI <- DTI$pasat[43:382]
x_DTI <- DTI$cca[43:382, ]
valid_indices <- which(!apply(x_DTI, 1, function(x) any(is.na(x))))
y_DTI <- y_DTI[valid_indices]
x_DTI <- x_DTI[valid_indices, ]
evaluate_and_plot(x_DTI, y_DTI, m_n1 = 20, n = 334, m_M0 = 100, nt = 93, nn = 93, file_name = "DTI.pdf")

