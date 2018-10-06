## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(qtl2pleio)
library(microbenchmark)

X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- gemma2::stagger_mats(X1, X1)
Sigma_inv <- diag(200)
Y <- runif(200)
calc_Bhat(X, Sigma_inv, Y)

microbenchmark(
  calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = Y), 
  rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = Sigma_inv)
)


