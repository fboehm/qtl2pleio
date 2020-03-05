library(qtl2pleio)
library(testthat)
context("testing calc_lrt_tib")

rep(rep(paste0('Marker', 1:2), times = 2), each = 2) -> marker1
rep(rep(paste0('Marker', 1:2), each = 2), each = 2) -> marker2
rep(paste0("Marker", 1:2), times = 4) -> marker3
runif(8, -1, 0) -> ll
tibble::tibble(marker1, marker2, ll) -> foo
bar <- tibble::tibble(marker1, marker2, marker3, ll)



test_that("calc_lrt_tib outputs a single number", {
  expect_true(is.numeric(calc_lrt_tib(foo)))
})

test_that("calc_lrt_tib works with 3-dim scan", {
  expect_true(is.numeric(calc_lrt_tib(bar)))
})

