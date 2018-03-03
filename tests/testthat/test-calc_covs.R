library(qtl2pleio)

context("calc_covs testing")

phe <- matrix(runif(200), ncol = 2)
sex <- matrix(rbinom(100, size = 1, prob = 1 / 2), ncol = 1)

calc_covs(phe, kinship = diag(100)) -> cc1
calc_covs(phe, kinship = diag(100), covariates = sex) -> cc2

phe3 <- matrix(runif(300), ncol = 3)


test_that("calc_covs returns different answers when covariates are used vs. not used", {
  expect_false(identical(cc1$Vg, cc2$Vg))
  expect_false(identical(cc1$Ve, cc2$Ve))
})


test_that("calc_covs accommodates d-variate phenotype for d more than 2", {
  expect_equal(calc_covs(phe3, diag(100))$Vg, calc_covs(phe3, diag(100))$Ve)
})
