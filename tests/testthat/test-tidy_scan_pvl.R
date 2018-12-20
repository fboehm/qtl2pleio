library(qtl2pleio)
library(testthat)
context("Testing tidy_scan_pvl function")

# set up dataframe

Var1 <- as.character(rep(paste0("m", 1:5), each = 5))
Var2 <- as.character(rep(paste0("m", 1:5), times = 5))
loglik <- - runif(25)
dat <- tibble::tibble(Var1, Var2, loglik)
pmap <- 1:5
names(pmap) <- as.character(paste0("m", 1:5))


test_that("tidy_scan_pvl outputs tibble with correct dimensions", {
  expect_equal(nrow(tidy_scan_pvl(dat, pmap)), as.integer(length(pmap) * 3))
})
