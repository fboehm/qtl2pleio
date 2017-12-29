library(qtl2pleio)
context("testing id_recombination")

test_that("logicals work as expected", {
  expect_true(id_recombination(c(0, 1), c(0, 1)),
  expect_false(id_recombination(c(1, 0), c(0, 1))))
})
