#' Create a bootstrap sample and perform multivariate QTL scan
#'
#' @param pp allele probabilities object for one chromosome only, like aprobs$`8`. Not a list
#' @param pleio_peak_index positive integer index indicating design matrix for simulation
#' @param phe n by d matrix of phenotypes
#' @param covariates n by n.cov matrix of covariates
#' @param kinship a kinship MATRIX, not a list
#' @param nboot_per_job number of bootstrap samples to call per invocation of function
#' @export
#' @return numerical vector of lrt statistics from one or more bootstrap samples
#'
boot_n <- function(pp, pleio_peak_index, phe, covariates = NULL, kinship, nboot_per_job = 1){
  X1 <- pp[ , , pleio_peak_index]
  if (!is.null(covariates)){cbind(X1, covariates) -> Xpre} else {X1 -> Xpre}
  ## remove subjects with missing values of phenotype
  is.na(phe[ , 1]) | is.na(phe[ , 2]) -> missing_indic
  phe_nona <- phe[!missing_indic, ]
  Xpre_nona <- Xpre[!missing_indic, ]
  k_nona <- k[!missing_indic, !missing_indic]
  ##
  gemma2::stagger_mats(Xpre_nona, Xpre_nona) -> X
  set.seed(proc_num)
  calc_covs(pheno = phe_nona, kinship = k_nona) -> cc_out
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
    scan_pvl(probs = pp[!missing_indic, , ], pheno = Ysim, covariates = sex[!missing_indic, drop = FALSE], kinship = k_nona, start_snp1 = s1, n_snp = nsnp) -> loglik
    # in above call, s1 & nsnp come from command line args
    calc_lrt_tib(loglik) -> lrt[i]
  }
  return(lrt)
}
