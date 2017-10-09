library(qtl2pleio)
context("testing scan_pvl")

# setup
## define probs
probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))
# define Y
Y_pre <- runif(200)
Y <- matrix(data = Y, nrow = 100)

test_that("pvl_scan returns a log-likelihood matrix", {
  expect_lte(scan_pvl(probs = probs, pheno = Y, kinship = diag(100), start_snp1 = 1, n_snp = 10)[1, 2], 0)
})
