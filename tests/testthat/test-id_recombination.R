library(qtl2pleio)
context("testing id_recombination in simple cases")


probs1 <- cbind(c(1, 1, 1, 0, 0, 0), c(0, 0, 0, 1, 1, 1))
probs2 <- cbind(c(1, 1, 1), c(0, 0, 0))


test_that("simple cases work as expected", {
  expect_true(id_recombination(probs1))
  expect_false(id_recombination(probs2))
})
