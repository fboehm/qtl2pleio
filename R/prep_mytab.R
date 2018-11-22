#' Prepare mytab object for use within scan_pvl R code
#'
#' @param d_size an integer, the number of traits
#' @param n_snp an integer, the number of markers
#' @export
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

