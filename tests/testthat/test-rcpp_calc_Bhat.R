library(qtl2pleio)
context("testing rcpp_calc_Bhat against calc_Bhat")

# setup
X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- gemma2::stagger_mats(X1, X1)
Sigma_inv <- diag(200)
Y <- runif(200)
# more complicated covariance matrix
Spre <- matrix(data = runif(200 * 200), nrow = 200, ncol = 200)
S <- Spre %*% t(Spre)
solve(S) -> S2_inv

test_that("rcpp_calc_Bhat and calc_Bhat give equal results", {
  expect_equal(
    rcpp_calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = Y),
    calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = Y)
  )
})

