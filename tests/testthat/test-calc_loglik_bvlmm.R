library(qtl2pleio)

context("testing calc_loglik_bvlmm")

Sigma <- diag(6)
B <- c(1, -1)
X <- gemma2::stagger_mats(rep(1, 3), rep(1, 3))
Y <- matrix(c(1:6), nrow = 3)
foo <- calc_loglik_bvlmm(X, B, Y, Sigma)
bar <- dnorm(as.vector(Y), mean = X %*% B, log = TRUE)
test_that("simple case - when independent - gives correct log likelihood value", {
  expect_equal(sum(bar), foo)
})
