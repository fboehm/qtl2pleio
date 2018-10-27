library(qtl2pleio)
library(testthat)
context("testing correct subsetting of inputs: phenotypes, allele probabilities, covariates, and kinship")

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

test_that("make_id2keep returns correct vector of subject ids both with and without covariates", {
  expect_length(make_id2keep(probs = probs, pheno = Y[1:10, ], addcovar = NULL, kinship = K1), 10)
  expect_length(make_id2keep(probs = probs, pheno = Y[1:10, ], addcovar = covariates, kinship = K1), 10)
})
