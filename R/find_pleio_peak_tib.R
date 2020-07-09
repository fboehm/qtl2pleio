#' Find the marker index corresponding to the peak of the pleiotropy trace in a tibble where the last column contains log likelihood values and the first d columns contain marker ids
#'
#' @param tib a (d+1) column tibble with first d columns containing marker ids and the last containing log likelihood values. Typically this is the output from `scan_pvl`.
#' @param start_snp positive integer, from the two-dimensional scan, that indicates where the scan started on the chromosome
#' @return positive integer indicating marker index for maximum value of log lik under pleiotropy
#' @export
#' @examples
#' marker1 <- rep(paste0('SNP', 1:3), times = 3)
#' marker2 <- rep(paste0('SNP', 1:3), each = 3)
#' loglik <- runif(9, -5, 0)
#' tibble::tibble(marker1, marker2, loglik) -> tib
#' find_pleio_peak_tib(tib, start_snp = 1)

find_pleio_peak_tib <- function(tib, start_snp) {
    nc <- ncol(tib)
    if (nc < 3){stop("tib must have at least 3 columns")}
    if (!is.numeric(unlist(tib[, nc]))){stop("last column of tib must be numeric")}
    smat <- as.matrix(tib[, -nc])
    # define an indicator vector for which rows to keep when making pleio_tib
    pleio_ind <- apply(FUN = function(x) {
        names(x) <- NULL
        identical(x, rep(x[1], nc - 1))
    }, X = smat, MARGIN = 1)
    pleio_tib <- tib[pleio_ind, ]
    ll <- unlist(pleio_tib[, nc])
    return(which.max(ll) + start_snp - 1)
}
