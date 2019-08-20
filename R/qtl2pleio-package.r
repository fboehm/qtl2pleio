## usethis namespace: start
#' @useDynLib qtl2pleio, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## usethis namespace: start
#' @import RcppEigen
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "lod", "marker_position"))

.onUnload <- function (libpath) {
  library.dynam.unload("qtl2pleio", libpath)
}
