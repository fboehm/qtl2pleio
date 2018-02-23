library(qtl2pleio)
context("testing scan_pvl")

# setup
## define probs
probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))
dimnames(probs)[[3]] <- paste0("Marker", 1:10)
# define Y
Y_pre <- runif(200)
Y <- matrix(data = Y_pre, nrow = 100)
covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
cov2 <- matrix(c(covariates[1:99], 10), nrow = 100, ncol = 1)
Y2 <- Y
Y2[1, 2] <- NA


scan_out <- scan_pvl(probs = probs, pheno = Y, kinship = diag(100), start_snp1 = 1, n_snp = 10)

test_that("pvl_scan returns a dataframe where the last column has numeric entries, all negative", {
  expect_true(identical(rep(TRUE, nrow(scan_out)), as.vector(scan_out[ , ncol(scan_out)] < 0)))
})

so_cov <- scan_pvl(probs = probs,
         pheno = Y, covariates = covariates,
         kinship = diag(100), start_snp1 = 1, n_snp = 10)

test_that("pvl_scan handles missing values in covariates appropriately", {
  expect_message(scan_pvl(probs = probs,
                          pheno = Y, covariates = covariates,
                          kinship = diag(100), start_snp1 = 1, n_snp = 10),
                 regexp = "1 subjects"
                 )
  expect_message(scan_pvl(probs = probs,
                          pheno = Y2, covariates = covariates,
                          kinship = diag(100), start_snp1 = 1, n_snp = 10),
                 regexp = "2 subjects")
  expect_equal(sum(!is.na(scan_out)), prod(dim(scan_out)))
  expect_equal(sum(!is.na(so_cov)), prod(dim(so_cov)))

})
