library(qtl2pleio)
library(testthat)
context("testing boot_pvl")

## define probs
probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))
s_id <- paste0('s', 1:100)
rownames(probs) <- s_id
colnames(probs) <- 'A'
dimnames(probs)[[3]] <- paste0('Marker', 1:10)
# define Y
Y_pre <- runif(200)
Y <- matrix(data = Y_pre, nrow = 100)
rownames(Y) <- s_id
colnames(Y) <- paste0('t', 1:2)
addcovar <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
rownames(addcovar) <- s_id
colnames(addcovar) <- 'c1'
kin <- diag(100)
rownames(kin) <- s_id
colnames(kin) <- s_id
Y2 <- Y
Y2[1, 2] <- NA
set.seed(2018-10-22)

test_that("output is vector of length nboot_per_job", {
  expect_length(boot_pvl(probs = probs, pheno = Y, kinship = kin,
           start_snp = 1, n_snp = 10, pleio_peak_index = 5,
           nboot_per_job = 1), 1)
  expect_length(boot_pvl(probs = probs, pheno = Y, kinship = kin,
                         start_snp = 1, n_snp = 10, pleio_peak_index = 5,
                         nboot_per_job = 4), 4)

})

test_that("output is numeric vector", {
  expect_true(is.numeric(boot_pvl(probs = probs, pheno = Y, kinship = kin,
                         start_snp = 1, n_snp = 10, pleio_peak_index = 5,
                         nboot_per_job = 1)))
  expect_true(is.numeric(boot_pvl(probs = probs, pheno = Y, kinship = kin,
                                  start_snp = 1, n_snp = 10, pleio_peak_index = 5,
                                  nboot_per_job = 4)))
})

