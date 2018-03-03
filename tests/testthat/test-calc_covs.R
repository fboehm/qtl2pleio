library(qtl2pleio)

context("calc_covs testing")

phe <- matrix(runif(200), ncol = 2)
sex <- matrix(rbinom(100, size = 1, prob = 1 / 2), ncol = 1)

calc_covs(phe, kinship = diag(100)) -> cc1
calc_covs(phe, kinship = diag(100), covariates = sex) -> cc2



test_that("calc_covs returns different answers when covariates are used vs. not used", {
  expect_false(identical(cc1$Vg, cc2$Vg))
  expect_false(identical(cc1$Ve, cc2$Ve))
})
