#' Check whether a vector, x, has all its entries equal to its first entry
#'
#' @param x a vector
#' @export
#' @return a logical indicating whether all vector entries are the same
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
#' An inputted matrix or 3-dimensional array is first subsetted - by rownames - to remove
#' those subjects who are not in `id2keep`. After that, the object's rows
#' are ordered to match the ordering of subject ids in the vector `id2keep`. This
#' (possibly reordered) object is returned.
#'
#' @param input a matrix of either phenotypes or covariates or array of allele probabilities
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
#' @examples
#' # define s_id
#' s_id <- paste0("s", 1:10)
#' # set up input matrix
#' foo <- matrix(data = rnorm(10 * 3), nrow = 10, ncol = 3)
#' rownames(foo) <- s_id
#' subset_input(input = foo, id2keep = s_id)
#' @return an object resulting from subsetting of `input`. Its rows are ordered per `id2keep`
#'
subset_input <- function(input, id2keep) {
    if (is.null(dim(input))){
        stop("input must have dim = 2 or 3")
    }
    if (length(dim(input)) == 2) {
        input <- input[rownames(input) %in% id2keep, , drop = FALSE]
        out <- input[match(rownames(input), id2keep), , drop = FALSE]
    }
   if (length(dim(input)) == 3) {
       input <- input[rownames(input) %in% id2keep, , , drop = FALSE]
       out <- input[match(rownames(input), id2keep), , , drop = FALSE]
    }
    return(out)
}

#' Subset a kinship matrix to include only those subjects present in all inputs
#'
#' Since a kinship matrix has subject ids in both rownames and colnames, so we need to
#' remove rows and columns according to names in `id2keep`. We first remove rows and
#' columns of subjects that are not in `id2keep`. We then order rows and columns of the
#' resulting matrix by the ordering in `id2keep`.
#'
#' @param kinship a kinship matrix
#' @param id2keep a character vector of subject ids to identify those subjects that are shared by all inputs
#' @export
subset_kinship <- function(kinship, id2keep) {
    if (is.null(dim(kinship))){
        stop("kinship matrix must have dim = 2")
    }
    kinship <- kinship[rownames(kinship) %in% id2keep, colnames(kinship) %in% id2keep]
    return(kinship[match(rownames(kinship), id2keep), match(colnames(kinship), id2keep)])
}

#' Identify shared subject ids among all inputs: covariates, allele probabilities array, kinship, and phenotypes
#'
#' We consider only those inputs that are not NULL. We then use `intersect` on pairs
#' of inputs' rownames to identify those subjects are shared among all non-NULL inputs.
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
#' We use `is.finite` from base R to identify those subjects that have one or
#' more missing values in `input_matrix`. We then return a character vector of subjects
#' that have no missingness in `input_matrix`.
#'
#' @param input_matrix phenotypes or covariates matrix
#' @export
#' @return character vector of subjects that have no missingness
check_missingness <- function(input_matrix){
    nr <- nrow(input_matrix)
    nc <- ncol(input_matrix)
    foo <- is.finite(input_matrix)
    bar <- matrix(data = foo, ncol = nc, nrow = nr)
    indic <- as.logical(apply(FUN = prod, X = bar, MARGIN = 1))
    return(rownames(input_matrix)[indic])
}


# drop linearly dependent columns
# if intercept=TRUE, add intercept before checking and then remove afterwards
# from Karl Broman's rqtl/qtl2 R package
#' @importFrom stats complete.cases
drop_depcols <-
    function(covar=NULL, add_intercept=FALSE, tol=1e-12)
    {
        if(is.null(covar)) return(covar)

        if(!is.matrix(covar)) covar <- as.matrix(covar)
        if(add_intercept) covar <- cbind(rep(1, nrow(covar)), covar)

        if(ncol(covar) <= 1) return(covar)

        X <- covar[complete.cases(covar),,drop=FALSE]

        # deal with NAs by omitting those rows before
        indep_cols <- sort(find_lin_indep_cols(X, tol))

        if(add_intercept) {

            target_ncol <- length(indep_cols)

            n_it <- 0
            while(!(1 %in% indep_cols)) {
                # don't want to omit the first column (the intercept)
                # need to work harder...
                #  - drop one column at a time other the intercept
                #  - when you find a column that doesn't reduce the target number of columns, drop it
                #  - check again if the intercept is being retained; if not, repeat

                for(i in seq_len(ncol(X))[-1]) { # loop over all but the first column (the intercept)
                    indep_cols <- find_lin_indep_cols(X[,-i,drop=FALSE], tol)
                    if(length(indep_cols) == target_ncol) { # this one definitely dependent; omit it
                        X <- X[,-i,drop=FALSE]
                        break
                    }
                }

                # new set of linearly independent columns
                indep_cols <- sort(find_lin_indep_cols(X, tol))

                n_it <- n_it+1 # count number of iterations and bail if it's large
                if(n_it > ncol(covar)+5) { # something is seriously messed up
                    stop("Not able to find set of linearly independent columns")
                }
            }

            # now drop the intercept
            indep_cols <- indep_cols[-1]
        }
        if(length(indep_cols)==0) return(NULL)

        covar[, indep_cols, drop=FALSE]
    }
