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

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "lod", "marker_position"))
