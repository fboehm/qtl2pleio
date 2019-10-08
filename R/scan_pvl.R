#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' `scan_pvl` calculates log likelihood for d-variate phenotype model fits. Inputted parameter `start_snp` indicates where in the `probs` object to start the scan.
#'
#' The function first discards individuals with one or more missing phenotypes or missing covariates.
#' It then infers variance components, Vg and Ve. Both Vg and Ve
#' are d by d covariance matrices. It uses an expectation maximization algorithm, as
#' implemented in the `gemma2` R package. `gemma2` R package is an R implementation of the
#' GEMMA algorithm for multivariate variance component estimation (Zhou & Stephens 2014 Nature methods).
#' Note that variance components are fitted on a model that uses the d-variate phenotype
#' but contains no genetic information. This model does, however,
#' use the specified covariates (after dropping dependent columns
#' in the covariates matrix).
#' These inferred covariance matrices, \eqn{\hat{Vg}} and \eqn{\hat{Ve}},
#' are then used in subsequent model fitting via
#' generalized least squares.
#' Generalized least squares model fitting is applied to every d-tuple of
#' markers within the specified genomic region for `scan_pvl`.
#' For a single d-tuple of markers, we fit the model:
#' \deqn{vec(Y) = Xvec(B) + vec(G) + vec(E)} where
#' \deqn{G \sim MN(0, K, \hat{Vg})} and \deqn{E \sim MN(0, I, \hat{Ve})} where \eqn{MN} denotes the matrix-variate
#' normal distribution with three parameters: mean matrix, covariance among rows, and
#' covariance among columns. \eqn{vec} denotes the vectorization operation, ie, stacking by columns.
#' \eqn{K} is a kinship matrix, typically calculated by leave-one-chromosome-out methods.
#' \eqn{Y} is the n by d phenotypes matrix. \eqn{X} is a block-diagonal nd by fd matrix consisting of
#' d blocks each of dimension n by f. Each n by f block (on the diagonal) contains a matrix of
#' founder allele probabilities for the n subjects at a single marker. The off-diagonal blocks
#' have only zero entries.
#' The log-likelihood is returned for each model. The outputted object is a tibble with
#' d + 1 columns. The first d columns specify the markers used in the corresponding model fit, while
#' the last column specifies the log-likelihood value at that d-tuple of markers.
#'
#' @param probs an array of founder allele probabilities for a single chromosome
#' @param pheno a matrix of phenotypes
#' @param kinship a kinship matrix for one chromosome
#' @param addcovar a matrix, n subjects by c additive covariates
#' @param start_snp index of where to start the scan within probs
#' @param n_snp the number of (consecutive) markers to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @param n_cores number of cores to use for calculations
#' @export
#' @importFrom stats var
#' @references Knott SA, Haley CS (2000) Multitrait
#' least squares for quantitative trait loci detection.
#' Genetics 156: 899–911.
#'
#' Jiang C, Zeng ZB (1995) Multiple trait analysis
#' of genetic mapping for quantitative trait loci.
#' Genetics 140: 1111-1127.
#'
#' Zhou X, Stephens M (2014) Efficient multivariate linear
#' mixed model algorithms for genome-wide association studies.
#' Nature methods 11:407-409.
#'
#' Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen S, Yandell BS, Churchill GA (2018)
#' R/qtl2: software for mapping quantitative trait loci with high-dimensional data and
#' multi-parent populations. Biorxiv https://www.biorxiv.org/content/early/2018/09/12/414748.
#'
#'
#' @examples
#' # read data
#' n <- 50
#' pheno <- matrix(rnorm(2 * n), ncol = 2)
#' rownames(pheno) <- paste0("s", 1:n)
#' colnames(pheno) <- paste0("tr", 1:2)
#' probs <- array(dim = c(n, 2, 5))
#' probs[ , 1, ] <- rbinom(n * 5, size = 1, prob = 0.2)
#' probs[ , 2, ] <- 1 - probs[ , 1, ]
#' rownames(probs) <- paste0("s", 1:n)
#' colnames(probs) <- LETTERS[1:2]
#' dimnames(probs)[[3]] <- paste0("m", 1:5)
#' scan_pvl(probs = probs, pheno = pheno, kinship = NULL,
#' start_snp = 1, n_snp = 5, n_cores = 1)
#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log10 likelihood

scan_pvl <- function(probs,
                     pheno,
                     kinship = NULL,
                     addcovar = NULL,
                     start_snp = 1,
                     n_snp,
                     max_iter = 1e+04,
                     max_prec = 1 / 1e+08,
                     n_cores = 1
                     )
    {
    if (is.null(probs)) stop("probs is NULL")
    if (is.null(pheno)) stop("pheno is NULL")
    stopifnot(!is.null(rownames(probs)),
              !is.null(colnames(probs)),
              !is.null(dimnames(probs)[[3]]),
              !is.null(rownames(pheno)),
              n_snp > 0,
              start_snp > 0,
              start_snp + n_snp - 1 <= dim(probs)[3]
        )
    # check additional conditions when addcovar is not NULL
    if (!is.null(addcovar)) {
        stopifnot(!is.null(rownames(addcovar)),
                  !is.null(colnames(addcovar))
                  )
    }
    d_size <- ncol(pheno)  # d_size is the number of univariate phenotypes
    # force things to be matrices
    if(!is.matrix(pheno)) {
        pheno <- as.matrix(pheno)
        if(!is.numeric(pheno)) stop("pheno is not numeric")
    }
    if(is.null(colnames(pheno))){ # force column names
        colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))}
    if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.numeric(addcovar)) stop("addcovar is not numeric")
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *one or more* phenotypes
    # need to consider presence or absence of different inputs: kinship, addcovar
    id2keep <- make_id2keep(probs = probs,
                            pheno = pheno,
                            addcovar = addcovar,
                            kinship = kinship
                            )
    # remove - from id2keep vector - subjects with a missing phenotype or covariate
    pheno <- subset_input(input = pheno, id2keep = id2keep)
    subjects_phe <- check_missingness(pheno)
    id2keep <- intersect(id2keep, subjects_phe)

    if (!is.null(addcovar)) {
        addcovar <- subset_input(input = addcovar, id2keep = id2keep)
        subjects_cov <- check_missingness(addcovar)
        id2keep <- intersect(id2keep, subjects_cov)
    }
    if (!is.null(addcovar)) {
        addcovar <- drop_depcols(addcovar)
    }
    # Send messages if there are two or fewer subjects
    if (length(id2keep) == 0){stop("no individuals common to all inputs")}
    if (length(id2keep) <= 2){
        stop(paste0("only ", length(id2keep),
                    " common individual(s): ",
                    paste(id2keep, collapse = ": ")))
        }
    # subset inputs to get all without missingness
    probs <- subset_input(input = probs, id2keep = id2keep)
    pheno <- subset_input(input = pheno, id2keep = id2keep)
    if (d_size != Matrix::rankMatrix(pheno)) stop("Phenotypes matrix is not full rank. Input only full-rank phenotypes matrices.")

    if (!is.null(kinship)) {
        kinship <- subset_kinship(kinship = kinship, id2keep = id2keep)
    }
    if (!is.null(addcovar)) {
        addcovar <- subset_input(input = addcovar, id2keep = id2keep)
    }
    if (!is.null(kinship)){
        # covariance matrix estimation
        message(paste0("starting covariance matrices estimation with data from ", length(id2keep), " subjects."))
        # first, run gemma2::MphEM(), by way of calc_covs(), to get Vg and Ve
        cc_out <- calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec, covariates = addcovar)
        Vg <- cc_out$Vg
        Ve <- cc_out$Ve
        message("covariance matrices estimation completed.")
        # define Sigma
        Sigma <- calc_Sigma(Vg, Ve, kinship)
    }
    if (is.null(kinship)){
        # get Sigma for Haley Knott regression without random effect
        Ve <- var(pheno) # get d by d covar matrix
        Sigma <- calc_Sigma(Vg = NULL, Ve = Ve, n_mouse = nrow(pheno))
    }

    # define Sigma_inv
    Sigma_inv <- solve(Sigma)
    # prepare table of marker indices for each call of scan_pvl
    mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
    # set up parallel analysis
    list_result <- parallel::mclapply(
                                X = as.data.frame(t(mytab)),
                                FUN = fit1_pvl,
                                addcovar = addcovar,
                                probs = probs,
                                inv_S = Sigma_inv,
                                S = Sigma,
                                start_snp = start_snp,
                                pheno = pheno,
                                mc.cores = n_cores
                                )
    mytab$loglik <- unlist(list_result)
    marker_id <- dimnames(probs)[[3]][start_snp:(start_snp + n_snp - 1)]
    mytab2 <- tibble::as_tibble(apply(FUN = function(x) marker_id[x], X = mytab[, -ncol(mytab)], MARGIN = 2))
    mytab2$log10lik <- mytab$loglik / log(10)
    return(mytab2)
}
