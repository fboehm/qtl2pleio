library(qtl2pleio)
context("testing rcpp_calc_Bhat2 against rcpp_calc_Bhat")

X <- matrix(data = rbinom(size = 1, n = 800, prob = 0.5), nrow = 100, ncol = 8)
Y <- matrix(data = rnorm(200), nrow = 100, ncol = 2)
Sigma <- diag(100)
Sigma_inv <- solve(Sigma)
Sigma2 <- Sigma + matrix(nrow = 100, ncol = 100, data = 0.2)
Sigma2_inv <- solve(Sigma2)

test_that("rcpp_calc_Bhat2 returns same values as rcpp_calc_Bhat", {
  expect_equal(
    rcpp_calc_Bhat(X = X / 2, Sigma = Sigma / 2, Y = Y / 2),
    rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = Sigma_inv)
  )
  expect_equal(
    rcpp_calc_Bhat(X = X, Sigma = Sigma2, Y = Y),
    rcpp_calc_Bhat2(X = X, Y = as.matrix(as.vector(Y[1:100])), Sigma_inv = Sigma2_inv)
  )
})

