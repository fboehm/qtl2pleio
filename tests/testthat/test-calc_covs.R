library(qtl2pleio)

context("calc_covs testing")

phe <- matrix(runif(200), ncol = 2)
sex <- matrix(rbinom(100, size = 1, prob = 1 / 2), ncol = 1)

calc_covs(phe, kinship = diag(100)) -> cc1
calc_covs(phe, kinship = diag(100), covariates = sex) -> cc2

phe3 <- matrix(runif(300), ncol = 3)
phe5 <- matrix(runif(500), ncol = 5)


test_that("calc_covs returns different answers when covariates are used vs. not used", {
  expect_false(identical(cc1$Vg, cc2$Vg))
  expect_false(identical(cc1$Ve, cc2$Ve))
})

test_that("calc_covs, when kinship is identity, has output that sums to give same result as cov()", {
  expect_equal(round(cc1$Vg + cc1$Ve, 4), round(cov(phe), 4), )
})

test_that("calc_covs accommodates d-variate phenotype for d more than 2 & Vg = Ve when kinship is identity", {
  expect_equal(calc_covs(phe3, diag(100))$Vg, calc_covs(phe3, diag(100))$Ve)
  expect_equal(calc_covs(phe5, diag(100))$Vg, calc_covs(phe5, diag(100))$Ve)
  expect_equal(calc_covs(phe5, diag(100), covariates = sex)$Vg, calc_covs(phe5, diag(100), covariates = sex)$Ve)
})

test_that("calc_covs outputs are full rank matrices", {
  expect_equal(qr(calc_covs(phe3, diag(100))$Vg)$rank, 3)
  expect_equal(qr(calc_covs(phe3, diag(100))$Ve)$rank, 3)
  expect_equal(qr(calc_covs(phe5, diag(100))$Ve)$rank, 5)
  expect_equal(qr(calc_covs(phe5, diag(100))$Vg)$rank, 5)
  expect_equal(qr(calc_covs(phe3, diag(100), covariates = sex)$Vg)$rank, 3)
  expect_equal(qr(calc_covs(phe3, diag(100), covariates = sex)$Ve)$rank, 3)
  expect_equal(qr(calc_covs(phe5, diag(100), covariates = sex)$Ve)$rank, 5)
  expect_equal(qr(calc_covs(phe5, diag(100), covariates = sex)$Vg)$rank, 5)
})

test_that("calc_covs outputs are symmetric", {
  expect_true(isSymmetric(calc_covs(phe3, diag(100))$Vg, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe3, diag(100))$Ve, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe5, diag(100))$Vg, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe5, diag(100))$Ve, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe3, diag(100), covariates = sex)$Vg, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe3, diag(100), covariates = sex)$Ve, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe5, diag(100), covariates = sex)$Vg, tol = 0.001))
  expect_true(isSymmetric(calc_covs(phe5, diag(100), covariates = sex)$Ve, tol = 0.001))
})
