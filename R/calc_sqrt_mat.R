#' Calculate matrix square root for a covariance matrix
#'
#' @param A covariance matrix
#' @export
#'
calc_sqrt_mat <- function(A){
  rcppeigen_sqrt(A)
}


.onUnload <- function (libpath) {
  library.dynam.unload("qtl2pleio", libpath)
}
