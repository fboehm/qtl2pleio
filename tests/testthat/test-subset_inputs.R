library(qtl2pleio)

context("ensure that subsetting of inputs works as intended")

x <- matrix(runif(100), nrow = 10, ncol = 10)
rownames(x) <- paste0("mouse", 1:10)

y <- x[10:1, ]
rownames(y) <- rownames(x)[10:1]
z <- x[10:1, 1, drop = FALSE]
rownames(z)[10] <- "mouse11"

id2keep <- make_id2keep(x, y)

id2keep2 <- make_id2keep()

test_that("subset_inputs arranges subject names in same order", {
  expect_identical(rownames(subset_input(input = x, id2keep = id2keep)),
                   rownames(subset_input(input = y, id2keep = id2keep))
                   )
  expect_false(identical(rownames(subset_input(input = x, id2keep = id2keep)),
                         rownames(y)
                         )
               )
  expect_identical(x, subset_input(input = y, id2keep = id2keep))
  })


