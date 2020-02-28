#' Extract founder allele effects at a single marker from output of qtl2::scan1coef
#'
#' @param marker_index an integer indicating where in the `map` object the peak position (or position of interest) is located
#' @param allele_effects_matrix output of `qtl2::scan1coef` for a single chromosome
#' @param map a map object for the chromosome of interest
#' @param columns which columns to choose within the `allele_effects_matrix`. Default is 1:8 to reflect 8 founder alleles of Diversity Outbred mice
#' @export
#' @return a vector of 8 founder allele effects at a single marker
#' @examples
#' # set up allele effects matrix
#' ae <- matrix(dat = rnorm(100 * 8), ncol = 8, nrow = 100)
#' ae[, 8] <- - rowSums(ae[, 1:7])
#' colnames(ae) <- LETTERS[1:8]
#' rownames(ae) <- paste0(1, "_", 1:100)
#' # set up map
#' map <- 1:100
#' names(map) <- rownames(ae)
#' # call get_effects
#' get_effects(marker_index = 15, allele_effects_matrix = ae, map = map)
#' @return a vector of founder allele effects at a single marker
get_effects <- function(marker_index, allele_effects_matrix, map, columns = 1:8){
  marker_name <- names(map[marker_index])
  return(allele_effects_matrix[rownames(allele_effects_matrix) == marker_name, columns])
}
