#' Fit a model for a specified d-tuple of markers
#'
#' `fit1_pvl` uses several functions in the package qtl2pleio to fit the
#' linear mixed effects model for a single d-tuple of markers.
#' Creation of `fit1_pvl` - from code that originally resided in `scan_pvl`, enabled
#' parallelization via the `parallel` R package.
#'
#' @param indices a vector of indices for extracting elements of `probs` array
#' @param start_snp an integer to specify the index of the marker where the scan - in call to scan_pvl - starts. This argument is needed because `mytab` has only relative indices (relative to the `start_snp` marker)
#' @param probs founder allele probabilities array
#' @param addcovar additive covariates matrix
#' @param inv_S inverse covariance matrix for the vectorized phenotype
#' @param S covariance matrix for the vectorized phenotype, ie, the inverse of inv_S. By making this a function input, we avoid inverting the matrix many many times.
#' @param pheno a n by d phenotypes matrix
#' @export
#' @return a number, the log-likelihood for the specified model
#' @examples
#' n <- 50
#' pheno <- matrix(rnorm(2 * n), ncol = 2)
#' Vg <- diag(2)
#' Ve <- diag(2)
#' Sigma <- calc_Sigma(Vg, Ve, diag(n))
#' Sigma_inv <- solve(Sigma)
#' probs <- array(dim = c(n, 2, 5))
#' probs[ , 1, ] <- rbinom(n * 5, size = 1, prob = 0.2)
#' probs[ , 2, ] <- 1 - probs[ , 1, ]
#' mytab <- prep_mytab(d_size = 2, n_snp = 5)
#' fit1_pvl(mytab[1, ], start_snp = 1,
#' probs = probs, addcovar = NULL, inv_S = Sigma_inv,
#' S = Sigma,
#' pheno = pheno
#' )

fit1_pvl <- function(indices,
                     start_snp,
                     probs,
                     addcovar,
                     inv_S,
                     S,
                     pheno
                     ){
  if (is.na(indices[length(indices)])) indices <- indices[- length(indices)]
  X_list <- prep_X_list(indices = indices,
                        start_snp = start_snp,
                        probs = probs,
                        covariates = addcovar
  )
  X <- gemma2::stagger_mats(X_list)
  Bhat <- rcpp_calc_Bhat2(X = X,
                          Sigma_inv = inv_S,
                          Y = as.vector(as.matrix(pheno))
                          )
  mymu <- as.vector(X %*% Bhat)
  out <- rcpp_log_dmvnorm2(inv_S = inv_S,
                           mu = mymu,
                           x = as.vector(as.matrix(pheno)),
                           S = S
                           )
  return(out)
}
