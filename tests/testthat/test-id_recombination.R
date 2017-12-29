library(qtl2pleio)
context("testing id_recombination")

test_that("logicals work as expected", {
  expect_true(id_recombination(),
  expect_false(probs_equal(c(1, 0), c(0, 1))))
})
