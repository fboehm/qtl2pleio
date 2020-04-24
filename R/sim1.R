#' Simulate a single data set consisting of n subjects and 2 phenotypes for each
#'
#' @param X design matrix (incorporating genotype probabilities from two loci), dn by df
#' @param B a matrix of allele effects, f rows by d columns
#' @param Vg a genetic covariance matrix, d by d
#' @param Ve an error covariance matrix, d by d
#' @param kinship a chromosome-specific kinship matrix, n by n
#' @export
#' @examples
#' n_mouse <- 20
#' geno <- rbinom(n = n_mouse, size = 1, prob = 1 / 2)
#' X <- gemma2::stagger_mats(geno, geno)
#' B <- matrix(c(1, 2), ncol = 2, nrow = 1)
#' Vg <- diag(2)
#' Ve <- diag(2)
#' sim1(X, B, Vg, Ve, kinship = diag(n_mouse))

sim1 <- function(X, B, Vg, Ve, kinship) {
    stopifnot(nrow(kinship) == ncol(kinship), nrow(Vg) == nrow(Ve), ncol(Vg) == ncol(Ve), nrow(Vg) ==
        ncol(Vg), nrow(Ve) == ncol(Ve))
    d_size <- ncol(Vg)
    Sigma <- calc_Sigma(Vg, Ve, kinship)
    n_mouse <- nrow(kinship)
    error_vec <- MASS::mvrnorm(n = 1, mu = rep(0, d_size * n_mouse), Sigma = Sigma)
    Y <- X %*% as.vector(B) + error_vec
    return(Y)
}
