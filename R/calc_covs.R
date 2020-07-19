#' Calculate Vg and Ve from d-variate phenotype and kinship
#'
#' @param pheno n by d matrix of phenotypes
#' @param kinship a kinship matrix, n by n
#' @param X1pre n by c design matrix. c = 1 to ignore genotypes
#' @param max_iter maximum number of EM iterations
#' @param max_prec maximum precision for stepwise increments in EM algorithm
#' @param covariates a n by n.cov matrix of numeric covariates
#' @return a list with 2 named components, Vg and Ve. Each is a d by d covariance matrix.
#' @examples
#' calc_covs(pheno = matrix(data = rnorm(100), nrow = 50, ncol = 2), kinship = diag(50))
#' @export
calc_covs <- function(pheno, kinship, X1pre = rep(1, nrow(kinship)), max_iter = 1e+06, max_prec = 1/1e+08,
    covariates = NULL) {
    e_out <- gemma2::eigen2(kinship)
    U <- e_out$vectors
    eval <- e_out$values
    n_mouse <- nrow(kinship)
    if (is.null(covariates)) {
        X1 <- t(X1pre) %*% U
    } else {
        X1 <- t(cbind(X1pre, covariates)) %*% U
    }
    Y <- t(pheno) %*% U
    d <- ncol(pheno)
    # run MphEM with only a design matrix that contains only the intercept term & covariates, if any(and
    # not any genotype info)
    foo <- gemma2::MphEM(max_iter = max_iter, max_prec = max_prec, X = X1, Y = Y, eval = eval, V_g = diag(d),
        V_e = diag(d))
    Vg <- foo[[length(foo)]][[2]]
    Ve <- foo[[length(foo)]][[3]]
    out <- list(Vg = Vg, Ve = Ve)
    return(out)
}
