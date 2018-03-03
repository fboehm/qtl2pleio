library(qtl2pleio)

context("calc_lrt testing")

foo <- runif(2)
ll1 <- diag(foo)
ll2 <- matrix(runif(4), nrow = 2, ncol = 2)

test_that("calc_lrt returns non-negative values", {
  expect_equal(calc_lrt(ll1), 0)
  expect_gte(calc_lrt(ll2), 0)
})

test_that("calc_lrt returns a number (vector of length 1)", {
  expect_equal(length(calc_lrt(ll2)), 1)
  expect_true(is.numeric(calc_lrt(ll2)))
})

test_that("calc_lrt gives an error when input contains a NA value", {
  expect_error(calc_lrt(diag(NA, 2)))
  expect_error(calc_lrt(diag(c(NA, 2))))

})
