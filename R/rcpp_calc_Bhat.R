#' Estimate allele effects matrix, B hat, with Rcpp
#'
#' @param X dn by df block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Sigma_inv dn by dn inverse covariance matrix, often composed as the inverse of \eqn{K \otimes V_g + I_n \otimes V_e}
#' @param Y dn by 1 matrix, ie, a column vector, of d phenotypes' measurements
#' @return a df by 1 matrix of GLS-estimated allele effects
#' @export
#'
rcpp_calc_Bhat <- function(X, Sigma_inv, Y){
  out <- solve(t(X) %*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv %*% as.vector(as.matrix(Y))
  return(out)
}
