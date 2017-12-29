library(qtl2pleio)
context("testing probs_equal")

test_that("logicals work as expected", {
  expect_true(probs_equal(c(0, 1), c(0, 1)),
  expect_false(probs_equal(c(1, 0), c(0, 1))))
})
