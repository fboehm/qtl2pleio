library(qtl2pleio)
library(testthat)
context("testing calc_lrt_tib")

rep(paste0('Marker', 1:3), times = 3) -> marker1
rep(paste0('Marker', 1:3), each = 3) -> marker2
runif(9, -1, 0) -> ll
tibble::tibble(marker1, marker2, ll) -> foo
bar <- tibble::tibble(marker1, marker2, ll, ll + 1)



test_that("calc_lrt_tib outputs a single number", {
  expect_true(is.numeric(calc_lrt_tib(foo)))
})

test_that("calc_lrt_tib outputs a single number", {
  expect_error(calc_lrt_tib(bar))
})

