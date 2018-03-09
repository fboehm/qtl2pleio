#' Calculate the LRT statistic
#'
#' @param loglik_mat a matrix of log likelihood values outputted by scan_pvl()
#'
#' @name calc_lrt-deprecated
#' @usage calc_lrt(loglik_mat)
#' @seealso \code{\link{qtl2pleio-deprecated}}
#' @keywords internal
NULL

#' @rdname qtl2pleio-deprecated
#' @section \code{calc_lrt}:
#' For \code{calc_lrt}, use \code{\link{qtl2pleio::calc_lrt_tib}}.
#'

#' @export

calc_lrt <- function(loglik_mat){
  .Deprecated(new = "calc_lrt_tib", package = "qtl2pleio")
  stopifnot(nrow(loglik_mat) == ncol(loglik_mat),
            sum(is.na(loglik_mat)) == 0,
            sum(!is.finite(loglik_mat)) == 0
            )
  out <- max(loglik_mat) - max(diag(loglik_mat))
  return(out)
}
