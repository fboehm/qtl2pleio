library(qtl2pleio)
context("testing rcpp_log_dmvnorm2")



test_that("rcpp_log_dmvnorm returns same values as mvtnorm::dmvnorm when log = TRUE", {
  expect_identical(
    as.numeric(rcpp_log_dmvnorm2(inv_S = diag(10), mu = rep(0, 10), x = rep(0, 10), S = diag(10))),
    mvtnorm::dmvnorm(x = rep(0, 10), mean = rep(0, 10), sigma = diag(10), log = TRUE)
  )
})

