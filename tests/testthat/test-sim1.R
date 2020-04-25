library(qtl2pleio)
context("Simulation of one phenotype vector")

X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X2 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
gemma2::stagger_mats(X1, X2) -> X
Sigma <- calc_Sigma(Vg = diag(1, 2), Ve = diag(1, 2), kinship = diag(100))


test_that("sim1 outputs numeric vector", {
  expect_true(is.numeric(sim1(X = X, B = c(1, 2), Sigma = Sigma)))
})

