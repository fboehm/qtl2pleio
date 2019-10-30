#' Prepare mytab object for use within scan_pvl R code
#'
#' @param d_size an integer, the number of traits
#' @param n_snp an integer, the number of markers
#' @export
#' @return a data.frame with d_size + 1 columns and (n_snp)^d_size rows. Last column is NA and named loglik.
#' @examples
#' prep_mytab(2, 10)
prep_mytab <- function(d_size, n_snp) {
    mytab <- expand.grid(rep(list(1:n_snp), d_size))
    mytab$loglik <- NA
    return(mytab)
}


#' Create a list of component X matrices for input to stagger_mats, to ultimately create design matrix
#'
#' @param indices a vector of integers
#' @param start_snp an integer denoting the index (within genotype probabilities array) where the scan should start
#' @param probs a three-dimensional array of genotype probabilities for a single chromosome
#' @param covariates a matrix of covariates
#' @export
#' @return a list of design matrices, ultimately useful when constructing the (multi-locus) design matrix
#' @examples
#' pp <- array(rbinom(n = 200, size = 1, prob = 0.5), dim = c(10, 2, 10))
#' prep_X_list(1:3, 1, probs = pp, covariates = NULL)
prep_X_list <- function(indices, start_snp, probs, covariates) {
    if (!is.null(covariates)) {
        out <- lapply(X = as.list(indices), FUN = function(x) {
            index <- x + start_snp - 1
            foo <- cbind(probs[, , index], covariates)
            return(foo)
        })
    } else {
        out <- lapply(X = as.list(indices), FUN = function(x) {
            index <- x + start_snp - 1
            foo <- probs[, , index]
            return(foo)
        })
    }
    return(out)
}

