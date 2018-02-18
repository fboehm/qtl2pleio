library(qtl2pleio)

context("testing arrange_by_rownames")

x <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(x) <- paste0("mouse", 1:10)

y <- x[10:1, ]
z <- x[10:1, 1, drop = FALSE]

test_that("check_dimnames works as expected", {
  expect_identical(rownames(arrange_by_rownames(y, x)), rownames(x))
  expect_false(check_dimnames(y, x))
  expect_equal(arrange_by_rownames(y, x), x)
  expect_equal(class(arrange_by_rownames(z, x)), "matrix")
  expect_equal(dim(arrange_by_rownames(z, x)), c(10, 1))
})
