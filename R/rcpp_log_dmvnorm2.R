#' Calculate log likelihood for a multivariate normal
#'
#' @param inv_S inverse covariance matrix
#' @param mu mean vector
#' @param x data vector
#' @param S covariance matrix, ie, the inverse of inv_S
#' @export

rcpp_log_dmvnorm2 <- function(inv_S, mu, x, S)
{
  n <- length(x)
  diag <- rcppeigen_get_diag(S)
  return((-n/2)*log(2*pi)-sum(log(diag))-(1/2)*(x-mu)%*%inv_S%*%(x-mu))
}
