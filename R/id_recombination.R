#' Determine if a single subject's genotype probabilities at two markers are equal.
#'
#' @param g1 a vector of genotype probabilities for a single subject at marker one
#' @param g2 a vector of genotype probabilities for the same subject at marker two
#' @export
probs_equal <- function(g1, g2){
  prod(g1 == g2) -> out
  return(as.logical(out))
}


#' Determine if all genotype probability vectors for markers in an interval are the same for one subject.
#'
#' @param index1 an integer denoting the left-hand marker.
#' @param index2 an integer denoting the right-hand marker.
#' @param probs a genotype probabilities 2-dimensional array for the appropriate chromosome and a single subject
#' @details Checks genotype probability vector (at each marker in interval) for equality with the left-hand marker's genotype probabilities vector.
#' @export
id_recombination <- function(index1, index2, probs){
  foo <- logical(length = index2 - index1)
  i <- 1
  for (index in (index1 + 1): index2){
    probs_equal(probs[index1, ], probs[index, ]) -> foo[i]
    i <- i + 1
  }
  prod(foo) -> out
  return(!as.logical(out))
}
