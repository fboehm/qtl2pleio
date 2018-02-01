library(qtl2pleio)
context("Check calc_Bhat")

# setup
X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- gemma2::stagger_mats(X1, X1)
Sigma_inv <- matrix(data = rep(1, 200 * 200), nrow = 200)
Y <- runif(200)


test_that("error is returned when inputted Sigma_inv is not positive definite", {
  expect_error(calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = Y))
})
test_that("dimensions of returned object are correct", {
  expect_equal(nrow(calc_Bhat(X = X, Sigma_inv = diag(200), Y = Y)), 2)
  expect_equal(ncol(calc_Bhat(X = X, Sigma_inv = diag(200), Y = Y)), 1)
})
