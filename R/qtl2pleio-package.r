#' qtl2pleio.
#'
#' Testing pleiotropy vs. separate QTL in multiparental populations
#'
#' @useDynLib qtl2pleio, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats profile
"_PACKAGE"

.onUnload <- function (libpath) {
  library.dynam.unload("qtl2pleio", libpath)
}

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "lod", "marker_position", "profile_lod", "loglik", "Var1", "profile", "pleio_max", "marker", "trait", "log10lik", "null_log10lik"))
