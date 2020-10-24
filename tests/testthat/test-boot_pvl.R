library(qtl2pleio)
library(testthat)
context("testing boot_pvl")

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

set.seed(2018-10-22)
nb <- 2
test_that("output is vector of length nboot", {
  expect_length(boot_pvl(probs = probs,
                         pheno = pheno,
                         kinship = kinship,
                         start_snp = 1,
                         n_snp = 10,
                         pleio_peak_index = 5,
                         nboot = nb, cores = 1),
                nb
                )

})
