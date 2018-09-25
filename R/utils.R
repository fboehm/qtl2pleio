#' Check whether a vector, x, has all its entries equal to its first entry
#'
#' @param x a vector
#' @export
#' @examples
#' x <- 1:5
#' check_identical(x)
#' y <- rep(1, 5)
#' check_identical(y)

check_identical <- function(x){
  stopifnot(is.vector(x))
  length(x) -> n
  return(identical(x, rep(x[1], n)))
}

#' Subset an input object - allele probabilities array or phenotypes matrix or covariates matrix. Kinship has its own subset function
#'
#' @param input a matrix of either phenotypes or covariates or array of allele probabilities
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
subset_input <- function(input, id2keep){
  if (length(dim(input)) == 2){
    out <- input[rownames(input) %in% id2keep, , drop = FALSE]
  }
  else {
    out <- input[rownames(input) %in% id2keep, , , drop = FALSE]
  }
  return(out)
}

#' Subset a kinship matrix to include only those subjects present in all inputs
#'
#' @param kinship a kinship matrix
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
subset_kinship <- function(kinship, id2keep){
  return(kinship[rownames(kinship) %in% id2keep, colnames(kinship) %in% id2keep])
}

#' Identify shared subject ids among all inputs: covariates, allele probabilities array, kinship, and phenotypes
#'
#' @param probs an allele probabilities array
#' @param pheno a phenotypes matrix
#' @param covar a covariates matrix
#' @param kinship a kinship matrix
#' @export
make_id2keep <- function(probs, pheno, covar, kinship){
  if (is.null(covar)){
    shared <- intersect(rownames(probs), rownames(pheno))
    id2keep <- intersect(shared, rownames(kinship))
  } else {
    shared1 <- intersect(rownames(probs), rownames(pheno))
    shared2 <- intersect(rownames(kinship), rownames(covar))
    id2keep <- intersect(shared1, shared2)
  }
  return(id2keep)
}
