#' qtl2pleio.
#'
#' @useDynLib qtl2pleio, .registration = TRUE
#' @name qtl2pleio
#' @docType package
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("qtl2pleio", libpath)
}
