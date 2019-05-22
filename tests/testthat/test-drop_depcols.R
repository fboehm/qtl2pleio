## The code here was originally prepared by Karl Broman
## https://github.com/rqtl/qtl2/blob/68883790388f996639f81b1bb2b7c244f5492164/tests/testthat/test-drop_depcols.R
## I edited the qtl2 test file to remove tests for functions that I don't use in qtl2pleio

library(qtl2pleio)

context("drop linearly dependent columns")

test_that("find_lin_indep_cols works", {

  set.seed(20151130)
  x <- cbind(1, sample(0:1, 200, repl=TRUE))

  expect_equal(sort(find_lin_indep_cols(x)), c(1,2))
  expect_equal(sort(find_lin_indep_cols(cbind(x, 1))), c(1,2))
  expect_equal(sort(find_lin_indep_cols(cbind(x, x[,1] + 0.5*x[,2]))), c(2,3))

  X <- matrix(rnorm(1000*10), ncol=10)
  expect_equal(sort(find_lin_indep_cols(X)), 1:10)
  expect_equal(length(find_lin_indep_cols(cbind(rowSums(X), X))), 10)

})

test_that("drop_depcols works", {

  set.seed(20151130)
  n.ind <- 100
  x <- cbind(1,
             sample(0:1, n.ind, replace=TRUE),
             sample(0:1, n.ind, replace=TRUE))

  expect_equal(drop_depcols(NULL), NULL)
  expect_equal(drop_depcols(x), x)
  for(i in 1:ncol(x))
    expect_equal(drop_depcols(x[,i]), x[,i,drop=FALSE])

  X <- cbind(rowSums(x), x)
  expect_equal(ncol(drop_depcols(X)), 3)

})


test_that("Fixed issue 22 re drop_depcols ", {

  # https://github.com/rqtl/qtl2/issues/22
  # ...case where indicators for all categories included

  # create example
  n_per <- 20
  n_group <- 4
  n <- n_per*n_group
  group <- rep(1:n_group, each=n_per)
  covar <- cbind(sex=rep(0:1,n_per*n_group),
                 diet=sample(100:200, n, replace=TRUE))
  for(i in 1:n_group)
    covar <- cbind(covar, (group==i)*1)
  colnames(covar)[-(1:2)] <- paste0("group", 1:n_group)

  # drop dependent columns
  covar_indep <- drop_depcols(covar, TRUE)

  # check that one of the groups got dropped
  expect_true("sex" %in% colnames(covar_indep))
  expect_true("diet" %in% colnames(covar_indep))
  expect_true(ncol(covar_indep) == ncol(covar)-1)

})
