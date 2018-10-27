library(qtl2pleio)
context("testing scan_pvl")

# setup
## define subject ids
s_id <- paste0("s", 101:200)
## define probs
probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))
dimnames(probs)[[3]] <- paste0("Marker", 1:10)
colnames(probs) <- "A"
rownames(probs) <- s_id
# define Y
Y_pre <- runif(200)
Y <- matrix(data = Y_pre, nrow = 100)
colnames(Y) <- c("y1", "y2")
rownames(Y) <- s_id
# define covariates
covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
colnames(covariates) <- "c1"
rownames(covariates) <- s_id
cov2 <- matrix(c(covariates[1:99], 10), nrow = 100, ncol = 1)
colnames(cov2) <- "c1"
rownames(cov2) <- s_id

Y2 <- Y
Y2[1, 2] <- NA
# define kinship
K1 <- diag(100)
rownames(K1) <- s_id
colnames(K1) <- s_id

## first scan_pvl call
scan_out <- scan_pvl(probs = probs,
                     pheno = Y,
                     kinship = K1,
                     start_snp = 1,
                     n_snp = 10
                     )

test_that("scan_pvl returns a dataframe where the last column has numeric entries, all negative", {
  expect_true(identical(rep(TRUE, nrow(scan_out)),
                        as.vector(scan_out[ , ncol(scan_out)] < 0)))
  expect_true(is.data.frame(scan_out))
})

so_cov <- scan_pvl(probs = probs,
                   pheno = Y,
                   addcovar = covariates,
                   kinship = K1,
                   start_snp = 1,
                   n_snp = 10
                   )

test_that("scan_pvl handles missing values in covariates appropriately", {
  expect_equal(sum(!is.na(scan_out)), prod(dim(scan_out)))
  expect_equal(sum(!is.na(so_cov)), prod(dim(so_cov)))
})

addcovar_rep <- matrix(rep(1, 200), ncol = 2)
rownames(addcovar_rep) <- s_id

test_that("scan_pvl handles linearly dependent columns in covariates by throwing error message", {
  expect_error(scan_pvl(probs = probs,
                        pheno = Y2,
                        addcovar = addcovar_rep,
                        kinship = K1,
                        start_snp = 1,
                        n_snp = 10
                        )
               )
}
)

