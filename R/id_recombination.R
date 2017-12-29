#' Determine if a single subject's lineage had a recombination (historically) between two markers
#'
#' @param g1 a vector of genotype probabilities for a single subject at marker one
#' @param g2 a vector of genotype probabilities for the same subject at marker two
#' @export
id_recombination <- function(g1, g2){
  prod(g1 == g2) -> out
  return(as.logical(out))
}
