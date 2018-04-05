#' Create a bootstrap sample and perform multivariate QTL scan
#'
#' @param pp allele probabilities object for one chromosome only, like aprobs$`8`. Not a list
#' @param pleio_peak_index positive integer index indicating design matrix for simulation
#' @param phe n by d matrix of phenotypes
#' @param covariates n by n.cov matrix of covariates
#' @param kinship a kinship matrix, not a list
#' @param nboot_per_job number of bootstrap samples to call per invocation of function
#' @param start_snp1 positive integer indicating index within pp for start of scan
#' @param n_snp number of markers to use in scan
#' @export
#' @return numerical vector of lrt statistics from one or more bootstrap samples
#'
boot_pvl <- function(pp, pleio_peak_index, phe, covariates = NULL, kinship, nboot_per_job = 1, start_snp1, n_snp){
  stopifnot(identical(nrow(kinship), nrow(pp)),
            identical(nrow(pp), nrow(phe)),
            check_dimnames(kinship, pp),
            check_dimnames(pp, phe)
            )
  if (!is.null(covariates)){
    stopifnot(check_dimnames(kinship, covariates),
              identical(nrow(pp), nrow(covariates))
              )
  }
  X1 <- pp[ , , pleio_peak_index]
  if (!is.null(covariates)){cbind(X1, covariates) -> Xpre} else {X1 -> Xpre}
  ## remove subjects with missing values of phenotype
  apply(FUN = function(x)identical(x, rep(TRUE, length(x))), X = !is.na(phe), MARGIN = 1) -> missing_indic
  phe_nona <- phe[!missing_indic, ]
  Xpre_nona <- Xpre[!missing_indic, ]
  k_nona <- k[!missing_indic, !missing_indic]
  ##
  gemma2::stagger_mats(Xpre_nona, Xpre_nona) -> X
  calc_covs(pheno = phe_nona, kinship = k_nona, covariates = covariates) -> cc_out
  (cc_out$Vg -> Vg)
  (cc_out$Ve -> Ve)
  # calculate Sigma
  calc_Sigma(Vg = Vg, Ve = Ve, K = k_nona) -> Sigma
  solve(Sigma) -> Sigma_inv
  # calc Bhat
  B <- calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = phe_nona)
  # Start loop
  lrt <- numeric()
  for (i in 1:nboot_per_job){
    sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = k_nona) -> foo
    matrix(foo, ncol = 2, byrow = FALSE) -> Ysim
    rownames(Ysim) <- rownames(phe_nona)
    colnames(Ysim) <- c("t1", "t2")
    if (!is.null(covariates)){
      scan_pvl(probs = pp[!missing_indic, , ], pheno = Ysim, covariates = covariates[!missing_indic, , drop = FALSE],
             kinship = k_nona, start_snp1 = s1, n_snp = nsnp) -> loglik
    } else {
      scan_pvl(probs = pp[!missing_indic, , ], pheno = Ysim,
               kinship = k_nona, start_snp1 = s1, n_snp = nsnp) -> loglik
    }
    # in above call, s1 & nsnp come from command line args
    calc_lrt_tib(loglik) -> lrt[i]
  }
  return(lrt)
}
