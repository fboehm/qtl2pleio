## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(qtl2pleio)
pheno <- matrix(rnorm(300), nrow = 100)
kinship <- diag(100)
calc_covs(pheno, kinship)

## ------------------------------------------------------------------------
library(limmbo2)
limmbo2(kinship = kinship, pheno = pheno, seed = 1)

## ------------------------------------------------------------------------
pheno <- matrix(rnorm(500), nrow = 100)
system.time(calc_covs(pheno, kinship))

## ------------------------------------------------------------------------
system.time(limmbo2(kinship = kinship, pheno = pheno, seed = 1))

## ------------------------------------------------------------------------
pheno <- matrix(rnorm(2500), nrow = 500)
kinship <- diag(500)
system.time(calc_covs(pheno, kinship) -> cc_out)
cc_out

## ------------------------------------------------------------------------
system.time(limmbo2(kinship = kinship, pheno = pheno, seed = 1) -> li_out)
li_out

## ------------------------------------------------------------------------
pheno <- matrix(rnorm(7500), nrow = 500)
system.time(calc_covs(pheno, kinship) -> cc_out)
cc_out

## ------------------------------------------------------------------------
system.time(limmbo2(kinship = kinship, pheno = pheno, S = 3, seed = 1) -> li_out)
li_out

