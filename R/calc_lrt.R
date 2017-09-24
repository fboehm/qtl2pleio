#' Calculate the LRT statistic
#'
#' @param loglik_mat a matrix of log likelihood values outputted by scan_pvl()
#' @export

calc_lrt <- function(loglik_mat){
  out <- max(loglik_mat) - max(diag(loglik_mat))
  return(out)
}
