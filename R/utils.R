#' Check whether a vector, x, has all its entries equal to its first entry
#'
#' @param x a vector
#' @returns a logical vector of length one
#' @export
#' @examples
#' x <- 1:5
#' check_identical(x)
#' y <- rep(1, 5)
#' check_identical(y)

check_identical <- function(x){
  stopifnot(is.vector(x))
  length(x) -> n
  return(all.equal(x, rep(x[1], n)))
}
