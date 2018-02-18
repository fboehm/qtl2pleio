#' Arrange rows of a matrix or data frame so that its rownames match those of another matrix or data frame
#'
#' @param df the object whose rows will be reordered
#' @param ref the object whose rows will not be reordered
#' @return df, but with rows in different order, so that the rownames match those of ref
#' @export
arrange_by_rownames <- function(df, ref){
  stopifnot(
    nrow(df) == nrow(ref)
  )
  df[match(x = dimnames(ref)[[1]], table = dimnames(df)[[1]]), , drop = FALSE] -> foo
  return(foo)
}

#' Verify that a single dimension's dimnames of two objects are equal
#'
#' @param x first object
#' @param y second object
#' @param xdim dimension to examine for first object
#' @param ydim dimension to examine for second object
#' @return logical indicator of equality of two objects' rownames
#' @export
check_dimnames <- function(x, y, xdim = 1, ydim = 1){
  identical(dimnames(x)[[xdim]], dimnames(y)[[ydim]])
}
