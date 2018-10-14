library(qtl2pleio)
context("testing rcpp_calc_Bhat2 and rcpp_calc_Bhat against calc_Bhat")

# setup
X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- gemma2::stagger_mats(X1, X1)
Sigma_inv <- diag(200)
Y <- runif(200)
# more complicated covariance matrix
Spre <- matrix(data = runif(200 * 200), nrow = 200, ncol = 200)
S <- Spre %*% t(Spre)
solve(S) -> S2_inv


test_that("rcpp_calc_Bhat2 returns same values as rcpp_calc_Bhat and calc_Bhat", {
  expect_equal(
    rcpp_calc_Bhat(X = X , Sigma_inv = Sigma_inv, Y = Y),
    rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = Sigma_inv)
  )
  expect_equal(
    rcpp_calc_Bhat(X = X , Sigma_inv = S2_inv, Y = Y),
    rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = S2_inv)
  )
  expect_equal(
    calc_Bhat(X = X , Sigma_inv = S2_inv, Y = Y),
    rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = S2_inv)
  )
})

