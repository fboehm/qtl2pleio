#' Calculate a LRT statistic from the (new) output from scan_pvl
#'
#' @param tib
#' @example
#' rep(paste0("Marker", 1:3), times = 3) -> marker1
#' rep(paste0("Marker", 1:3), each = 3) -> marker2
#' runif(9, -1, 0) -> ll
#' tibble::tibble(marker1, marker2, ll) -> scan_out
#' calc_lrt_tib(scan_out)
#' @export
#' @return a number, the LRT statistic
calc_lrt_tib <- function(tib){
  ncol(tib) -> nc
  # define an indicator vector for which rows to keep when making pleio_tib
  pleio_ind <- apply(FUN = function(x) identical(x, rep(x[1], nc - 1)), X = tib[ , - nc], MARGIN = 1)
  pleio_tib <- tib[pleio_ind, ]
  return(max(tib$ll) - max(pleio_tib$ll))
}
