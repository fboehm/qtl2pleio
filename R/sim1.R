#' Simulate a single data set consisting of n subjects and 2 phenotypes for each
#'
#' @param X design matrix (incorporating genotype probabilities from two loci), 2n by 2a
#' @param B a matrix of allele effects, a rows by 2 columns
#' @param Vg a genetic covariance matrix, 2 by 2
#' @param Ve an error covariance matrix, 2 by 2
#' @param kinship a chromosome-specific kinship matrix, n by n
#' @export

sim1 <- function(X, B, Vg, Ve, kinship){
  stopifnot(nrow(kinship) == ncol(kinship), nrow(Vg) == nrow(Ve),
            ncol(Vg) == ncol(Ve),
            nrow(Vg) == ncol(Vg),
            nrow(Ve) == ncol(Ve)
            )
  nrow(X) / 2 -> n_mouse
  In <- diag(n_mouse)
  Sigma <- kinship %x% Vg + In %x% Ve
  error_vec <- MASS::mvrnorm(n = 1, mu = rep(0, 2 * n_mouse), Sigma = Sigma)
  Y <- X %*% as.matrix(as.vector(B)) + error_vec
  return(Y)
}
