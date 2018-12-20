library(qtl2pleio)
library(testthat)
context("Testing add_pmap function")

# set up tibble

marker1 <- as.character(rep(paste0("m", 1:5), each = 5))
marker2 <- as.character(rep(paste0("m", 1:5), times = 5))
loglik <- - runif(25)
dat <- tibble::tibble(marker1, marker2, loglik)
pmap <- 1:5
names(pmap) <- as.character(paste0("m", 1:5))

test_that("add_pmap outputs tibble with correct dimensions", {
  expect_equal(nrow(add_pmap(dat, pmap)), nrow(dat))
  expect_equal(ncol(add_pmap(dat, pmap)), ncol(dat) + 2)
  expect_true(tibble::is_tibble(add_pmap(dat, pmap)))
})
