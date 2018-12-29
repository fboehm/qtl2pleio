#' Calculate a likelihood ratio test statistic from the output of scan_pvl()
#'
#' @param tib outputted tibble from scan_pvl
#' @examples
#' rep(paste0('Marker', 1:3), times = 3) -> marker1
#' rep(paste0('Marker', 1:3), each = 3) -> marker2
#' runif(9, -1, 0) -> ll
#' tibble::tibble(marker1, marker2, ll) -> scan_out
#' calc_lrt_tib(scan_out)
#' @export
#' @return a number, the (log) likelihood ratio test statistic
calc_lrt_tib <- function(tib) {
    nc <- ncol(tib)
    if (nc != 3) stop("tib must have exactly 3 columns")
    smat <- as.matrix(tib[, -nc])
    # define an indicator vector for which rows to keep when making pleio_tib
    pleio_ind <- apply(FUN = function(x) {
        names(x) <- NULL
        identical(x, rep(x[1], nc - 1))
    }, X = smat, MARGIN = 1)
    pleio_tib <- tib[pleio_ind, ]
    return(max(tib[, nc]) - max(pleio_tib[, nc]))
}
