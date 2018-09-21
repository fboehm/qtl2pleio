#' Estimate allele effects matrix, B hat, with Rcpp
#'
#' @param X dn by df block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Sigma dn by dn covariance matrix, often composed as \eqn{K \otimes V_g + I_n \otimes V_e}
#' @param Y dn by 1 matrix, ie, a column vector, of d phenotypes' measurements
#' @return a df by 1 matrix of GLS-estimated allele effects
#' @export
#'
rcpp_calc_Bhat <- function(X, Sigma, Y){
  calc_invsqrt_mat(Sigma) -> invsqrt_Sigma
  Wmat <- invsqrt_Sigma %*% X
  Zmat <- invsqrt_Sigma %*% Y
  out <- rcpp_ols(X = Wmat, Y = Zmat)
  return(out)
}

#' Estimate allele effects matrix, B hat, with Rcpp
#'
#' @param X dn by df block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Y dn by 1 matrix, ie, a column vector, of d phenotypes' measurements
#' @param Sigma_inv dn by dn inverse covariance matrix, often composed as inverse of \eqn{K \otimes V_g + I_n \otimes V_g}
#' @return a df by 1 matrix of GLS-estimated allele effects
#' @export
#'
rcpp_calc_Bhat2 <- function(X, Y, Sigma_inv){
  out <- rcpp_gls(X = X, Y = Y, Sigma_inv = Sigma_inv)
  return(out)
}

