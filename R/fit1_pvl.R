#' Fit a model for a specified d-tuple of markers
#'
#' `fit1_pvl` uses several functions in the package qtl2pleio to fit the
#' linear mixed effects model for a single d-tuple of markers.
#' Creation of `fit1_pvl` - from code that originally resided in `scan_pvl`, enabled
#' parallelization via the `parallel` R package.
#'
#' @param mytab a `mytab` object, with d+ 1 columns, typically the output of `prep_mytab`
#' @param rownum an integer to indicate which row of `mytab` to use
#' @param probs founder allele probabilities matrix
#' @param addcovar additive covariates matrix
#' @param Sigma_inv inverse covariance matrix for the vectorized phenotype
#' @param pheno a n by d phenotypes matrix
#' @export
#' @return a number, the log-likelihood for the specified model
#' @examples
#' # read data
#' iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' # insert pseudomarkers into map
#' map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
#' # calculate genotype probabilities
#' probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' # leave-one-chromosome-out kinship matrices
#' kinship <- qtl2::calc_kinship(probs, "loco")
#' # get founder allele probabilites
#' probs <- qtl2::genoprob_to_alleleprob(probs)
#' cc_out <- calc_covs(pheno, kinship$`1`, max_iter = max_iter, max_prec = max_prec, covariates = addcovar)
#' Vg <- cc_out$Vg
#' Ve <- cc_out$Ve
#' # define Sigma
#' Sigma <- calc_Sigma(Vg, Ve, kinship)
#' # define Sigma_inv
#' Sigma_inv <- solve(Sigma)
#' # prepare table of marker indices for each call of scan_pvl
#' mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
#' # set up parallel analysis
#' fit1_pvl(mytab, 1, start_snp = 1,
#' probs = probs$`1`, addcovar = NULL, Sigma_inv = Sigma_inv,
#' pheno = pheno
#' )

fit1_pvl <- function(mytab, rownum,
                     start_snp,
                     probs,
                     addcovar,
                     Sigma_inv,
                     pheno
                     ){
  indices <- unlist(mytab[rownum, ])
  X_list <- prep_X_list(indices = indices[-length(indices)],
                        start_snp = start_snp,
                        probs = probs,
                        covariates = addcovar
  )
  X <- gemma2::stagger_mats(X_list)
  Bhat <- rcpp_calc_Bhat2(X = X,
                          Sigma_inv = Sigma_inv,
                          Y = as.vector(as.matrix(pheno))
                          )
  mymu <- as.vector(X %*% Bhat)
  out <- rcpp_log_dmvnorm2(inv_S = Sigma_inv,
                           mu = mymu,
                           x = as.vector(as.matrix(pheno)),
                           S = Sigma
                           )
  return(out)
}
