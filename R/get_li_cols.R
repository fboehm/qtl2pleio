#' Extract a subset of covariate columns that are linearly independent
#'
#'
#'
#' We have the option of adding or not adding an intercept before doing the QR decomposition,
#' as indicated by the input `add_intercept`.
#' We will typically use `add_intercept = TRUE` because, in practice,
#' `find_li_cols` is called within `scan_pvl`, the function that
#' performs the two-dimensional scan. Within `scan_pvl`, we form a design matrix
#' that consists of a founder allele probabilities matrix AND the covariates.
#' To ensure that this matrix is invertible, we need to consider the possibility
#' appending the covariates to the founder allele probabilities matrix leads
#' to a matrix that is not full rank, possibly due to the fact that
#' the founder allele probabilities add to 1 and that a subset of the covariate columns
#' may sum to 1 (or a multiple of 1).
#'
#' @references
#' https://stackoverflow.com/questions/14943422/r-independent-columns-in-matrix
#' @examples
#' c1 <- rep(c(0, 1), each = 4)
#' c2 <- rep(c(0, 1), times = 4)
#' c3 <- rep(c(0, 0, 1, 1), times = 2)
#' addcovar <- cbind(c1, c2, c3)

get_li_cols <- function(addcovar, add_intercept = TRUE){
  # check if each column has only one value repeated, eg, (2, 2, 2, ..., 2)
  id_logical <- apply(FUN = function(x)identical(x, rep(x[1], length(x))),
                      X = addcovar,
                      MARGIN = 2
                      )
  addcovar <- addcovar[ , !id_logical] #remove column(s) that have a single value only
  rr <- qr(addcovar)$rank # rank of inputted matrix, after possibly removing column(s) with all entries a single value
  if (add_intercept){
    addcovar <- cbind(1, addcovar)
  }
  q <- qr(addcovar)
  out <- addcovar[ , q$pivot[seq(q$rank)], drop = FALSE]
  if (add_intercept){ # remove the column of all 1's before returning the result
    out <- out[ , -1, drop = FALSE]
  }
  return(out)
}
