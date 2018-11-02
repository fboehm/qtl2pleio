library(qtl2pleio)
library(testthat)

context("testing get_li_cols")

c1 <- rep(c(0, 1), each = 4)
c2 <- rep(c(0, 1), times = 4)
c3 <- rep(c(0, 0, 1, 1), times = 2)
ac1 <- cbind(c2, c3, c1)
ac2 <- cbind(ac1, 1)

test_that("get_li_cols returns a matrix of same rank as inputted matrix,
          with `add_intercept = TRUE`", {
            expect_equal(qr(get_li_cols(ac1))$rank, qr(ac1)$rank)
            expect_equal(qr(get_li_cols(ac2))$rank, qr(ac2)$rank - 1)
          }
          )
test_that("get_li_cols returns a matrix of same rank as inputted matrix,
          with `add_intercept = FALSE`", {
            expect_equal(qr(get_li_cols(ac1,
                                        add_intercept = FALSE))$rank,
                         qr(ac1)$rank
            )
            expect_equal(qr(get_li_cols(ac2,
                                        add_intercept = FALSE))$rank,
                         qr(ac2)$rank - 1)
            # rank - 1 because the first column is all 1's and we need to discard it typically
          }
)
