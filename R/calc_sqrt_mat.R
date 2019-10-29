#' Calculate matrix square root for a covariance matrix
#'
#' @param A covariance matrix
#'
calc_sqrt_mat <- function(A) {
    rcppeigen_sqrt(A)
}


#' Calculate matrix inverse square root for a covariance matrix
#'
#' @param A covariance matrix
#'
calc_invsqrt_mat <- function(A) {
    rcppeigen_invsqrt(A)
}



