library(qtl2pleio)

context("testing check_dimnames")

x <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(x) <- paste0("mouse", 1:10)
colnames(x) <- rownames(x)
y <- x
z <- x
rownames(z) <- paste0("frog", 1:10)
colnames(z) <- rownames(z)

test_that("check_dimnames works as expected", {
  expect_true(check_dimnames(x, y))
  expect_false(check_dimnames(x, z))
})
