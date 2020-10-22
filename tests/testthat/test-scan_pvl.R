library(qtl2pleio)
library(testthat)
context("testing scan_pvl")

# setup
# read data
is_inst <- function(pkg) {
  nzchar(system.file(package = pkg))
}
qtl2_indic <- is_inst("qtl2")

if (qtl2_indic) {
  iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
  # insert pseudomarkers into map
  map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
  # calculate genotype probabilities
  probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
  # grab phenotypes and covariates; ensure that covariates have names attribute
  pheno <- iron$pheno
  # leave-one-chromosome-out kinship matrices
  kinship <- qtl2::calc_kinship(probs, "loco")$`1`
  # get founder allele probabilites
  probs <- qtl2::genoprob_to_alleleprob(probs)$`1`
  ac <- matrix(as.numeric(iron$covar$sex == "m", ncol = 1))
  colnames(ac) <- "sex"
  rownames(ac) <- rownames(probs)

  ## first scan_pvl call
  scan_out <- scan_pvl(probs = probs,
                       pheno = pheno,
                       kinship = kinship,
                       start_snp = 1,
                       n_snp = 10, cores = 1
                       )

  test_that("scan_pvl returns a dataframe where the last column has numeric entries, all negative", {
    expect_true(identical(rep(TRUE, nrow(scan_out)),
                          as.vector(scan_out[ , ncol(scan_out)] < 0)))
    expect_true(is.data.frame(scan_out))
  })

  so_cov <- scan_pvl(probs = probs,
                     pheno = pheno,
                     addcovar = ac,
                     kinship = kinship,
                     start_snp = 1,
                     n_snp = 10, cores = 1
                     )



  foo <- scan_pvl(probs = probs,
           pheno = pheno,
           kinship = NULL,
           start_snp = 1,
           n_snp = 10, cores = 1
  )

  test_that("scan_pvl handles missing values in covariates appropriately", {
    expect_equal(sum(!is.na(scan_out)), prod(dim(scan_out)))
    expect_equal(sum(!is.na(so_cov)), prod(dim(so_cov)))
  })

  test_that("scan_pvl handles kinship = NULL", {
    expect_true(is.data.frame(foo))
  })



  test_that("scan_pvl throws error when phenotypes matrix is not full rank",{
    pheno_singular <- cbind(pheno[ , 1], 2 * pheno[ , 1])
    expect_error(scan_pvl(probs = probs,
                          pheno = pheno_singular,
                          kinship = NULL,
                          start_snp = 1,
                          n_snp = 10, cores = 1
    ),
    regexp = "Phenotypes matrix is not full rank")
  })
}
