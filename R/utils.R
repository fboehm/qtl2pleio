#' Check whether a vector, x, has all its entries equal to its first entry
#'
#' @param x a vector
#' @export
#' @examples
#' x <- 1:5
#' check_identical(x)
#' y <- rep(1, 5)
#' check_identical(y)

check_identical <- function(x) {
    stopifnot(is.vector(x))
    n <- length(x)
    return(identical(x, rep(x[1], n)))
}

#' Subset an input object - allele probabilities array or phenotypes matrix or covariates matrix. Kinship has its own subset function
#'
#' @param input a matrix of either phenotypes or covariates or array of allele probabilities
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
subset_input <- function(input, id2keep) {
    if (is.null(dim(input))){
        stop("input must have dim = 2 or 3")
    }
    if (length(dim(input)) == 2) {
        out <- input[match(rownames(input), id2keep), , drop = FALSE]
    }
   if (length(dim(input)) == 3) {
        out <- input[match(rownames(input), id2keep), , , drop = FALSE]
    }
    return(out)
}

#' Subset a kinship matrix to include only those subjects present in all inputs
#'
#' @param kinship a kinship matrix
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
subset_kinship <- function(kinship, id2keep) {
    if (is.null(dim(kinship))){
        stop("kinship matrix must have dim = 2")
    }
    return(kinship[match(rownames(kinship), id2keep), match(colnames(kinship), id2keep)])
}

#' Identify shared subject ids among all inputs: covariates, allele probabilities array, kinship, and phenotypes
#'
#' @param probs an allele probabilities array
#' @param pheno a phenotypes matrix
#' @param addcovar a covariates matrix
#' @param kinship a kinship matrix
#' @export
#' @return a character vector of subject IDs common to all (non-null) inputs
make_id2keep <- function(probs,
                         pheno,
                         addcovar = NULL,
                         kinship = NULL) {
    if (is.null(addcovar) & !is.null(kinship)) {
        shared <- intersect(rownames(probs), rownames(pheno))
        id2keep <- intersect(shared, rownames(kinship))
    }
    if (!is.null(addcovar) & !is.null(kinship)){
        shared1 <- intersect(rownames(probs), rownames(pheno))
        shared2 <- intersect(rownames(kinship), rownames(addcovar))
        id2keep <- intersect(shared1, shared2)
    }
    if (is.null(addcovar) & is.null(kinship)) {
        id2keep <- intersect(rownames(probs), rownames(pheno))
    }
    if (!is.null(addcovar) & is.null(kinship)){
        shared1 <- intersect(rownames(probs), rownames(pheno))
        id2keep <- intersect(shared1, rownames(addcovar))
    }
    return(id2keep)
}

#' Check for missingness in phenotypes or covariates
#'
#' @param input_matrix phenotypes or covariates matrix
#' @export
#' @return character vector of subjects that have no missingness
check_missingness <- function(input_matrix){
    nr <- nrow(input_matrix)
    nc <- ncol(input_matrix)
    foo <- input_matrix[is.finite(input_matrix)]
    bar <- matrix(data = foo, ncol = nc, nrow = nr)
    indic <- as.logical(apply(FUN = prod, X = bar, MARGIN = 1))
    return(rownames(input_matrix)[indic])
}


