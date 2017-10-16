#' Calculate estimated allele effects, B matrix
#'
#' @param X 2n by 2f block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Sigma_inv 2n by 2n inverse covariance matrix, often composed as the inverse of \eqn{K \otimes V_g + I_n \otimes V_e}
#' @param Y 2n by 1 matrix, ie, a column vector, of 2 phenotypes' measurements
#' @export
calc_Bhat <- function(X, Sigma_inv, Y){
  out <- solve(t(X) %*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv %*% as.vector(Y)
  return(out)
}
