#' Calculate the phenotypes covariance matrix Sigma
#'
#' @param Vg 2 by 2 genetic covariance matrix for the two phenotypes
#' @param Ve 2 by 2 error covariance matrix for the two phenotypes
#' @param K n by n kinship matrix
#' @export
calc_Sigma <- function(Vg, Ve, K){
  n_mouse <- nrow(K)
  out <- Vg %x% K  + Ve %x% diag(n_mouse)
  return(out)
}
