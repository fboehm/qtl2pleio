library(qtl2pleio)
context("testing prep_mytab output's structure")


test_that("prep_mytab produces output with correct dimensions", {
  expect_equal(ncol(prep_mytab(d_size = 3, n_snp = 4)), 4)
  expect_equal(ncol(prep_mytab(d_size = 4, n_snp = 4)), 5)
  expect_equal(nrow(prep_mytab(d_size = 3, n_snp = 4)), 4 ^ 3)
  expect_equal(nrow(prep_mytab(d_size = 4, n_snp = 4)), 4 ^ 4)
})


test_that("verify that last column is all NAs", {
  expect_identical(as.vector(prep_mytab(d_size = 2, n_snp = 4)[, 3]),
                   rep(NA, 4 ^ 2))
})
