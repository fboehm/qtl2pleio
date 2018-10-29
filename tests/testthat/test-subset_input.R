library(qtl2pleio)
library(testthat)
context("testing correct subsetting - based on common ids - of inputs (with proper
        ordering of subject ids): phenotypes, allele probabilities,
        covariates, and kinship")

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


# tests

test_that("subset_input, with covariates, when not NULL,
          returns a matrix with the correct number of rows and correct
          subject ids in row names", {
  expect_equal(nrow(subset_input(input = covariates,
                                id2keep = paste0("s", 101:110))),
               10)
  expect_equal(ncol(subset_input(input = covariates,
                                           id2keep = paste0("s", 101:110))),
                         ncol(covariates))

  expect_identical(rownames(subset_input(input = covariates,
                                    id2keep = paste0("s", 101:110))),
                   paste0("s", 101:110))
  expect_identical(rownames(subset_input(input = covariates,
                                         id2keep = paste0("s", 110:101))),
                   paste0("s", 110:101))

})

test_that("subset_input, with pheno, returns
          a matrix with the correct number of
          rows and correct subject ids in row names", {
  expect_equal(nrow(subset_input(input = Y,
                                id2keep = paste0("s", 101:110)
  )), 10)
  expect_equal(ncol(subset_input(input = Y,
                                 id2keep = paste0("s", 101:110)
            )), ncol(Y))
  expect_identical(rownames(subset_input(input = Y,
                                         id2keep = paste0("s", 101:110)
                                         )), paste0("s", 101:110))
  expect_true(is.matrix(subset_input(input = Y,
                                     id2keep = paste0("s", 101:110))))
})

test_that("subset_input, with allele probabilities array,
          returns a matrix with the correct number of rows
          and correct subject ids in row names and that
          dimension = 3 is preserved", {
  expect_equal(nrow(subset_input(input = probs,
                                id2keep = paste0("s", 101:110)
                                )), 10)
  expect_identical(rownames(subset_input(input = probs,
                                         id2keep = paste0("s", 101:110)
                                         )), paste0("s", 101:110))
  expect_equal(length(dim(subset_input(input = probs,
                                         id2keep = paste0("s", 101:110)
                                       ))), 3)
})

test_that("subset_kinship returns a matrix with the correct number of
          rows & columns and correct subject ids in row names and column names", {
  expect_equal(dim(subset_kinship(kinship = K1,
                                  id2keep = paste0("s", 101:110)))[1],
               10)
  expect_equal(dim(subset_kinship(kinship = K1,
                                  id2keep = paste0("s", 101:110))
                   )[2],
               10)
  expect_identical(rownames(subset_kinship(kinship = K1,
                                         id2keep = paste0("s", 101:110))),
                   paste0("s", 101:110))
  expect_identical(colnames(subset_kinship(kinship = K1,
                                         id2keep = paste0("s", 101:110)
  )), paste0("s", 101:110))
  expect_identical(colnames(subset_kinship(kinship = K1,
                                           id2keep = paste0("s", 110:101)
  )), paste0("s", 110:101))
  expect_identical(rownames(subset_kinship(kinship = K1,
                                           id2keep = paste0("s", 110:101)
  )), paste0("s", 110:101))

})

x <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(x) <- paste0("mouse", 1:10)

y <- x[10:1, ]
rownames(y) <- rownames(x)[10:1]
z <- x[10:1, 1, drop = FALSE]
rownames(z)[10] <- "mouse11"

id2keep <- make_id2keep(probs = x, pheno = y)



test_that("subset_inputs arranges subject names in order specified in id2keep", {
  expect_identical(rownames(subset_input(input = x, id2keep = id2keep)),
                   rownames(subset_input(input = y, id2keep = id2keep))
  )
  expect_false(identical(rownames(subset_input(input = x, id2keep = id2keep)),
                         rownames(y)
  )
  )
  expect_identical(x, subset_input(input = y, id2keep = id2keep))
})





