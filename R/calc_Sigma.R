#' Calculate the phenotypes covariance matrix Sigma
#'
#' @param Vg d by d genetic covariance matrix for the d phenotypes
#' @param Ve d by d error covariance matrix for the d phenotypes
#' @param K n by n kinship matrix
#' @export
calc_Sigma <- function(Vg, Ve, K) {
    n_mouse <- nrow(K)
    out <- Vg %x% K + Ve %x% diag(n_mouse)
    return(out)
}
