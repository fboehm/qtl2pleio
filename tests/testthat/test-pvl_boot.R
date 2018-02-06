library(qtl2pleio)

context("testing parametric bootstrap code")



probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))

X <- as.matrix(probs[ , , 5])
nboot <- 5

test_that("pvl_boot returns a numeric vector of length nboot", {
  expect_length(pvl_boot(X = X, B = c(1, 2), Vg_initial = diag(2), Ve_initial = diag(2), kinship = diag(100), probs = probs, start_snp = 1, n_snp = 10, nboot = nboot), nboot)
  expect_true(is.numeric(pvl_boot(X = X, B = c(1, 2), Vg_initial = diag(2), Ve_initial = diag(2), kinship = diag(100), probs = probs, start_snp = 1, n_snp = 10, nboot = nboot)))
  expect_gte(pvl_boot(X = X, B = c(1, 2), Vg_initial = diag(2), Ve_initial = diag(2), kinship = diag(100), probs = probs, start_snp = 1, n_snp = 10, nboot = nboot)[1], 0)
})
