#' Calculate the phenotypes covariance matrix Sigma
#'
#' @param Vg d by d genetic covariance matrix for the d phenotypes
#' @param Ve d by d error covariance matrix for the d phenotypes
#' @param kinship optional n by n kinship matrix. if NULL, Vg is not used.
#' @param n_mouse number of subjects
#' @export
#' @return dn by dn covariance matrix
calc_Sigma <- function(Vg, Ve, kinship = NULL, n_mouse = nrow(kinship)) {
    if (!is.null(kinship)){
      out <- Vg %x% kinship + Ve %x% diag(n_mouse)
    }
    if (is.null(kinship)){
      out <- Ve %x% diag(n_mouse) # n_mouse must be specified when kinship = NULL
    }
    return(out)
}
