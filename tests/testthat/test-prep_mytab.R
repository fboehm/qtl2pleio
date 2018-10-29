library(qtl2pleio)
context("testing prep_mytab")


test_that("prep_mytab produces output with correct dimensions", {
  expect_equal(ncol(prep_mytab(d_size = 3, n_snp = 4)), 4)
  expect_equal(ncol(prep_mytab(d_size = 4, n_snp = 4)), 5)
  expect_equal(nrow(prep_mytab(d_size = 3, n_snp = 4)), 4 ^ 3)
  expect_equal(nrow(prep_mytab(d_size = 4, n_snp = 4)), 4 ^ 4)
})
