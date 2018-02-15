library(qtl2pleio)

context("testing check_dimnames")

x <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(x) <- paste0("mouse", 1:10)

y <- x[10:1, ]

test_that("check_dimnames works as expected", {
  expect_identical(rownames(arrange_by_rownames(y, x)), rownames(x))
  expect_false(check_dimnames(y, x))
  expect_equal(arrange_by_rownames(y, x), x)
})
