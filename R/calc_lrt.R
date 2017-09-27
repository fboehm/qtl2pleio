#' Calculate the LRT statistic
#'
#' @param loglik_mat a matrix of log likelihood values outputted by scan_pvl()
#' @export

calc_lrt <- function(loglik_mat){
  stopifnot(nrow(loglik_mat) == ncol(loglik_mat),
            sum(is.na(loglik_mat)) == 0,
            sum(!is.finite(loglik_mat)) == 0
            )
  out <- max(loglik_mat) - max(diag(loglik_mat))
  return(out)
}
