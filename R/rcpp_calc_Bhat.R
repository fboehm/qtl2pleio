#' Estimate allele effects matrix, B hat, with Rcpp functions
#'
#' @param X dn by df block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Sigma_inv dn by dn inverse covariance matrix, where its inverse, ie, Sigma, is often composed as \eqn{K \otimes V_g + I_n \otimes V_e}
#' @param Y dn by 1 matrix, ie, a column vector, of d phenotypes' measurements
#' @return a df by 1 matrix of GLS-estimated allele effects
#' @examples
#' X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
#' X <- gemma2::stagger_mats(X1, X1)
#' Sigma_inv <- diag(200)
#' Y <- runif(200)
#' rcpp_calc_Bhat(X = X, Sigma_inv = Sigma_inv, Y = Y)
#' @export
#'
rcpp_calc_Bhat <- function(X, Sigma_inv, Y) {
    # calc_invsqrt_mat(Sigma) -> invsqrt_Sigma
    invsqrt_Sigma <- calc_sqrt_mat(Sigma_inv)
    Wmat <- invsqrt_Sigma %*% X
    Zmat <- invsqrt_Sigma %*% Y
    out <- rcpp_ols(X = Wmat, Y = Zmat)
    return(out)
}

#' Estimate allele effects matrix, B hat, with Rcpp functions
#'
#' @param X dn by df block-diagonal design matrix that incorporates genetic info for two markers. Note that we can use the same marker data twice.
#' @param Y dn by 1 matrix, ie, a column vector, of d phenotypes' measurements
#' @param Sigma_inv dn by dn inverse covariance matrix, often composed as inverse of \eqn{K \otimes V_g + I_n \otimes V_g}
#' @return a df by 1 matrix of GLS-estimated allele effects
#' @examples
#' X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
#' X <- gemma2::stagger_mats(X1, X1)
#' Sigma_inv <- diag(200)
#' Y <- runif(200)
#' rcpp_calc_Bhat2(X = X, Y = Y, Sigma_inv = Sigma_inv)
#' @export
#'
rcpp_calc_Bhat2 <- function(X, Y, Sigma_inv) {
    out <- rcpp_gls(X = X, Y = Y, Sigma_inv = Sigma_inv)
    return(out)
}

