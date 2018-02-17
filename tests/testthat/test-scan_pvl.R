library(qtl2pleio)
context("testing scan_pvl")

# setup
## define probs
probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))
# define Y
Y_pre <- runif(200)
Y <- matrix(data = Y_pre, nrow = 100)
covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
cov2 <- matrix(c(covariates[1:99], 10), nrow = 100, ncol = 1)

scan_out <- scan_pvl(probs = probs, pheno = Y, kinship = diag(100), start_snp1 = 1, n_snp = 10)

test_that("pvl_scan returns a log-likelihood matrix with correctly named rows and cols", {
  expect_lte(scan_out[1, 2], 0)
  expect_equal(class(scan_out), "matrix")
  expect_equal(nrow(scan_out), ncol(scan_out))
  expect_identical(rownames(scan_out), colnames(scan_out))
})

so_cov <- scan_pvl(probs = probs,
         pheno = Y, covariates = covariates,
         kinship = diag(100), start_snp1 = 1, n_snp = 10)

test_that("pvl_scan handles missing values in covariates appropriately", {
  expect_message(scan_pvl(probs = probs,
                          pheno = Y, covariates = covariates,
                          kinship = diag(100), start_snp1 = 1, n_snp = 10), regexp = "1 subjects"
                 )
  expect_equal(sum(!is.na(scan_out)), prod(dim(scan_out)))
  expect_equal(sum(!is.na(so_cov)), prod(dim(so_cov)))

})
