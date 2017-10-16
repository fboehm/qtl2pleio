#' Calculate estimated allele effects, B matrix
#'
#' @param X 2n by 2f block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Sigma 2n by 2n covariance matrix, often composed as \eqn{K \otimes V_g + I_n \otimes V_e}
#' @param Y 2n by 1 matrix, ie, a column vector, of 2 phenotypes' measurements
#' @export
calc_Bhat <- function(X, Sigma, Y){
  stopifnot(det(Sigma) > 0, nrow(X) == nrow(Y))
  Sig_inv <- solve(Sigma)
  return(solve(t(X) %*% Sig_inv %*% X) %*% t(X) %*% Sig_inv %*% as.vector(Y))
}
