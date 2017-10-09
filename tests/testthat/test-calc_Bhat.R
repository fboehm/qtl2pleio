library(qtl2pleio)
context("Check inputs to calc_Bhat")

# setup
X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- pleiotropy::stagger_mats(X1, X1)
Sigma <- matrix(data = rep(1, 200 * 200), nrow = 200)
Y <- runif(200)


test_that("error is returned when inputted Sigma is not positive definite", {
  expect_error(calc_Bhat(X = X, Sigma = Sigma, Y = Y))
})
