library(qtl2pleio)
context("testing rcpp_calc_Bhat")

X <- matrix(data = rbinom(size = 1, n = 800, prob = 0.5), nrow = 100, ncol = 8)
Y <- matrix(data = rnorm(200), nrow = 100, ncol = 2)
Sigma <- diag(100)

test_that("rcpp_calc_Bhat returns same values as GLS with calculations in R", {
  expect_equal(
    rcpp_calc_Bhat(X = X, Sigma = Sigma, Y = Y),
    solve(t(X) %*% solve(Sigma) %*% X) %*% t(X) %*% solve(Sigma) %*% Y
  )
})

