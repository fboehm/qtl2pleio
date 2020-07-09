library(qtl2pleio)
library(testthat)
context("Testing find_pleio_peak_tib")

marker1 <- rep(paste0('SNP', 1:3), times = 3)
marker2 <- rep(paste0('SNP', 1:3), each = 3)
set.seed(2018-12-25)
loglik <- runif(9, -5, 0)
tibble::tibble(marker1, marker2, loglik) -> tib
tibble::tibble(marker1, marker2, 1:9) -> tib2
tibble::tibble(marker1, marker2, - loglik) -> tib3

test_that("find_pleio_peak_tib returns the index that has the greatest value of loglik",
          {
            expect_equal(which.max(loglik[c(1, 5, 9)]),
                         as.integer(find_pleio_peak_tib(tib, start_snp = 1)))
          })

test_that("Input checks for find_pleio_peak_tib work as expected",
          {
            expect_error(find_pleio_peak_tib(tib2[ , 1:2], start_snp = 1))
            })
