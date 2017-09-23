#' Calculate likelihood ratio test statistics for simulated data sets for a single marker
#'
#' @param X a genotypes probability matrix for a single marker
#' @param B allele effects matrix, often calculated via calc_Bhat
#' @param Vg_initial estimated Vg covariance matrix from true data
#' @param Ve_initial estimated Ve covariance matrix from true data
#' @param kinship a chromosome-specific kinship matrix
#' @param nboot number of bootstrap samples
#' @export

pvl_boot <- function(X, B, Vg_initial, Ve_initial, kinship, nboot = 1000){
  boot_stat <- numeric(length = nboot)
  n_mouse <- nrow(X)
  gemma2::eigen2(kinship) -> e_out
  e_out$values -> eval
  e_out$vectors -> U
  for (i in 1:nboot){
    # simulate phenotype data
    Y <- sim1(X = X, B = B, Vg = Vg_initial, Ve = Ve_initial, kinship = kinship)
    # estimate Vg & Ve for the simulated data
    gemma2::MphEM(eval = eval, X = t(X) %*% U, Y = t(Y) %*% U, V_g = diag(2), V_e = diag(2)) -> mph_out
    mph_out[[length(mph_out)]][[2]] -> Vg
    mph_out[[length(mph_out)]][[3]] -> Ve
    Sigma <- kinship %x% Vg + diag(n_mouse) %x% Ve
    Bvec <- calc_Bhat(pleiotropy::stagger_mats(X, X), Sigma, Y)
    Bmat <- matrix(data = Bvec, byrow = FALSE, ncol = 2)
    boot_stat[i] <- calc_loglik_bvlmm(X = X, Y = Y, Bhat = B, Sigma = Sigma)
  }
  return(boot_stat)
}
