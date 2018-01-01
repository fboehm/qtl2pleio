#' Determine if a single subject's genotype probabilities at two markers are equal.
#'
#' @param g1 a vector of genotype probabilities for a single subject at marker one
#' @param g2 a vector of genotype probabilities for the same subject at marker two
#' @export
probs_equal <- function(g1, g2){
  rg1 <- round(g1, digits = 3)
  rg2 <- round(g2, digits = 3)
  prod(rg1 == rg2) -> out
  return(as.logical(out))
}


#' Determine if all genotype probability vectors for markers in an interval are the same for one subject.
#'
#' @param probs a genotype probabilities 2-dimensional array for the appropriate chromosome and a single subject
#' @details Checks genotype probability vector (at each marker in probs) for equality with the left-hand marker's genotype probabilities vector.
#' @export
id_recombination <- function(probs){
  stopifnot(nrow(probs) > 1,
            ncol(probs) > 1
            )
  foo <- logical(length = nrow(probs) - 1)
  for (i in 2:nrow(probs)){
    probs_equal(probs[1, ], probs[i, ]) -> foo[i - 1]
  }
  prod(foo) -> out
  return(!as.logical(out))
}
