#' Convert `scan_multi_oneqtl` output of `qtl2::scan1` output
#'
#' We convert output of `scan_multi_oneqtl` into format outputted by `qtl2::scan1`.
#'
#' @param sm_output tibble output from scan_multi_oneqtl for one chromosome only
#' @param trait_name character vector (of length one) specifying the trait names
#' @return object of class `scan1`
#' @export
#' @examples
#' # read data
#' iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' \dontshow{iron <- iron[,7]}
#'
#' # insert pseudomarkers into map
#' map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
#'
#' # calculate genotype probabilities
#' probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
#'
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' covar <- match(iron$covar$sex, c("f", "m")) # make numeric
#' names(covar) <- rownames(iron$covar)
#' Xcovar <- qtl2::get_x_covar(iron)
#'
#' aprobs <- qtl2::genoprob_to_alleleprob(probs)
#' sm_out <- scan_multi_oneqtl(probs = aprobs, pheno = pheno)
#' sm_to_s1 <- convert_to_scan1_output(sm_out[[1]], trait_name = "tr1and2")
#'
#' # 95% Bayes credible interval for QTL on chr 7, first phenotype
#' qtl2::bayes_int(sm_to_s1, map)
#'
convert_to_scan1_output <- function(sm_output, trait_name){
  lods <- sm_output %>%
    dplyr::mutate(lod = log10lik - null_log10lik) %>%
    dplyr::select(Var1, lod)
  mat <- as.matrix(lods$lod)
  rownames(mat) <- lods$Var1
  colnames(mat) <- trait_name
  class(mat) <- c("scan1", "matrix")
  return(mat)
}
