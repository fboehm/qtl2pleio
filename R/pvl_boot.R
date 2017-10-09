#' Calculate likelihood ratio test statistics for simulated data sets for a single marker
#'
#' @param X a genotypes probability matrix for a single marker, ie, the marker that has the peak LOD under pleiotropy
#' @param B allele effects matrix, often calculated via calc_Bhat
#' @param Vg_initial estimated Vg covariance matrix from true data
#' @param Ve_initial estimated Ve covariance matrix from true data
#' @param kinship a chromosome-specific kinship matrix
#' @param probs a single chromosome's genotype probabilities
#' @param start_snp index for starting the scan
#' @param n_snp number of SNPs to include in the scan
#' @param nboot number of bootstrap samples
#' @export

pvl_boot <- function(X, B, Vg_initial, Ve_initial, kinship, probs, start_snp, n_snp, nboot = 1000){
  boot_stat <- numeric(length = nboot)
  n_mouse <- nrow(X)
  X_2 <- pleiotropy::stagger_mats(X, X)
  for (i in 1:nboot){
    # simulate phenotype data
    Y <- sim1(X = X_2, B = B, Vg = Vg_initial, Ve = Ve_initial, kinship = kinship)
    loglik_mat <- scan_pvl(probs = probs, pheno = Y, kinship = kinship, start_snp1 = start_snp, n_snp = n_snp)
    boot_stat[i] <- calc_lrt(loglik_mat)
  }
  return(boot_stat)
}
