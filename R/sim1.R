#' Simulate a single data set consisting of n subjects and 2 phenotypes for each
#'
#' @param X genotype probabilities matrix for a single marker
#' @param B a matrix of allele effects
#' @param Vg a genetic covariance matrix
#' @param Ve an error covariance matrix
#' @param kinship a chromosome-specific kinship matrix
#' @export

sim1 <- function(X, B, Vg, Ve, kinship){
  nrow(X) -> n_mouse
  In <- diag(n_mouse)
  Sigma <- kinship %x% Vg + In %x% Ve
  error_vec <- MASS::mvrnorm(n = 1, mu = rep(0, 2 * n_mouse), Sigma = Sigma)
  error_mat <- matrix(data = error_vec, byrow = FALSE, ncol = 2)
  Y <- X %*% B + error_mat
  return(Y)
}
