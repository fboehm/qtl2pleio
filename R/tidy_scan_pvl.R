#' Tidy the matrix of log likelihood values for further analysis & plotting
#'
#' @param loglik_mat a (square) matrix of log likelihood values
#' @export

tidy_scan_pvl <- function(loglik_mat){
  # Assumes that we have rownames and columns names assigned to loglik_mat
  marker1 <- rep(rownames(loglik_mat), times = ncol(loglik_mat))
  marker2 <- rep(colnames(loglik_mat), each = nrow(loglik_mat))
  ll <- as.vector(loglik_mat)
  return(tibble::tibble(marker1, marker2, ll))
}
