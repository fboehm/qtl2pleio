#' Simulate a single multivariate data set consisting of n subjects and d phenotypes for each
#'
#' @param X design matrix (incorporating genotype probabilities from two loci), dn by df
#' @param B a matrix of allele effects, f rows by d columns
#' @param Sigma dn by dn covariance matrix
#' @export
#' @return a vector of length dn. The first n entries are for trait 1, the second n for trait 2, etc.
#' @examples
#' n_mouse <- 20
#' geno <- rbinom(n = n_mouse, size = 1, prob = 1 / 2)
#' X <- gemma2::stagger_mats(geno, geno)
#' B <- matrix(c(1, 2), ncol = 2, nrow = 1)
#' sim1(X, B, Sigma = diag(2 * n_mouse))

sim1 <- function(X, B, Sigma) {
    error_vec <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    Y <- X %*% as.vector(B) + error_vec
    return(Y)
}
