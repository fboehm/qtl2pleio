#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' `scan_pvl` calculates log likelihood for d-variate phenotype model fits. Inputted parameter `start_snp` indicates where in the `probs` object to start the scan.
#'
#' The function first discards individuals with missing phenotypes or missing addcovar.
#' It also discards addcovar that have the same value for all remaining subjects.
#' It then uses one of two methods to infer variance components, Vg and Ve. Both Vg and Ve
#' are d by d covariance matrices.
#' The first method, the default, uses an expectation maximization algorithm, as
#' implemented in the `gemma2` R package. `gemma2` R package is an R implementation of the
#' GEMMA algorithm (Zhou & Stephens 2012 Nature Genetics).
#' The alternative method for inferring variance components, which is useful for more than two
#' phenotypes, ie, when d > 2, uses the LiMMBo method to infer variance components (Meyer et al. 2018).
#' Specifically, with input `use_limmbo2 = TRUE`, `scan_pvl` calls the function
#' `limmbo2::limmbo2` to infer d by d variance components.
#' Note that variance components, when inferred via
#' `gemma2`, are fitted on a model that uses the d-variate phenotype
#' but contains no genetic information. This model does, however, use the specified covariates.
#' Current implementation of `limmbo2` doesn't use covariates.
#' These inferred covariance matrices, \eqn{\hat{Vg}} and \eqn{\hat{Ve}}, are then used in subsequent model fitting via
#' generalized least squares. Generalized least squares model fitting is applied to every d-tuple of
#' markers within the specified genomic region for `scan_pvl`.
#' A progress bar appears *after* inferring Vg and Ve and records progress through the
#' d-tuples of markers that constitute the d-variate scan.
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
#'
#' @param probs an array of founder allele probabilities for a single chromosome
#' @param pheno a matrix of phenotypes
#' @param kinship a kinship matrix for one chromosome
#' @param addcovar a matrix, n subjects by n.cov, of additive covariates, where each column is one numeric additive covariate
#' @param start_snp index of where to start the scan within probs
#' @param n_snp the number of (consecutive) markers to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @export
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
#'
#'
#' @examples
#' ## define probs
#'probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
#'probs <- array(data = probs_pre, dim = c(100, 1, 10))
#'s_id <- paste0('s', 1:100)
#'rownames(probs) <- s_id
#'colnames(probs) <- 'A'
#'dimnames(probs)[[3]] <- paste0('Marker', 1:10)
#'# define Y
#'Y_pre <- runif(200)
#'Y <- matrix(data = Y_pre, nrow = 100)
#'rownames(Y) <- s_id
#'colnames(Y) <- paste0('t', 1:2)
#'covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
#'rownames(covariates) <- s_id
#'colnames(covariates) <- 'c1'
#'kin <- diag(100)
#'rownames(kin) <- s_id
#'colnames(kin) <- s_id
#'Y2 <- Y
#'Y2[1, 2] <- NA
#'scan_pvl(probs = probs, pheno = Y, kinship = kin,
#'         start_snp = 1, n_snp = 10)
#'scan_pvl(probs = probs, pheno = Y2, kinship = kin,
#'         start_snp = 1, n_snp = 10)
#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log likelihood

scan_pvl <- function(probs,
                     pheno,
                     kinship = NULL,
                     addcovar = NULL,
                     start_snp = 1,
                     n_snp,
                     max_iter = 1e+06,
                     max_prec = 1 / 1e+08
                     )
    {
    if(is.null(probs)) stop("probs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    stopifnot(!is.null(rownames(probs)),
              !is.null(colnames(probs)),
              !is.null(dimnames(probs)[[3]]),
              !is.null(rownames(pheno)),
              !is.null(colnames(pheno)),
              #!is.null(rownames(kinship)),
              #!is.null(colnames(kinship)),
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
    # check that the objects have rownames
    qtl2::check4names(pheno, kinship, probs)
    # force things to be matrices
    if(!is.matrix(pheno)) {
        pheno <- as.matrix(pheno)
        if(!is.numeric(pheno)) stop("pheno is not numeric")
    }
    if(is.null(colnames(pheno))) # force column names
        colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
    if(!is.null(addcovar)) {
        if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if(!is.numeric(addcovar)) stop("addcovar is not numeric")
    }

    # find individuals in common across all arguments
    # and drop individuals with missing covariates or missing *one or more* phenotypes
    ind2keep <- qtl2::get_common_ids(probs, addcovar, complete.cases=TRUE)
    ind2keep <- qtl2::get_common_ids(ind2keep, rownames(pheno)[rowSums(!is.finite(pheno)) > 0])
    if(length(ind2keep) <= 2) {
        if(length(ind2keep) == 0)
            stop("No individuals in common.")
        else
            stop("Only ", length(ind2keep), " individuals in common: ",
                 paste(ind2keep, collapse = ":"))
    }
    # make sure addcovar is full rank when we add an intercept
    addcovar <- qtl2:::drop_depcols(addcovar, TRUE, tol)



    # remove mice with missing values of phenotype or missing value(s) in addcovar
    missing_indic <- matrix(!apply(FUN = is.finite,
                                   X = pheno,
                                   MARGIN = 1
                                   ),
                            nrow = nrow(pheno),
                            ncol = ncol(pheno),
                            byrow = TRUE
                            )
    missing2 <- apply(FUN = function(x) {
        identical(as.logical(x),
                  rep(FALSE, ncol(pheno))
                  )
        },
        MARGIN = 1,
        X = missing_indic
        )
    if (!is.null(addcovar)) {
        addcovar <- subset_input(input = addcovar,
                                   id2keep = intersect(rownames(pheno),
                                                       rownames(addcovar))
                                   )
        pheno <- subset_input(input = pheno,
                              id2keep = intersect(rownames(pheno),
                                                  rownames(addcovar))
                              )
        miss_cov <- matrix(!apply(FUN = is.finite,
                                  X = addcovar,
                                  MARGIN = 1),
                           nrow = nrow(addcovar),
                           ncol = ncol(addcovar),
                           byrow = TRUE
                           )
        miss_cov2 <- apply(FUN = function(x){
            identical(as.logical(x), rep(FALSE, ncol(addcovar)))
            },
            MARGIN = 1,
            X = miss_cov
            )
        missing2 <- missing2 & miss_cov2
    }
    if (sum(!missing2) > 0) {
        message(paste0(sum(!missing2), " subjects dropped due to missing values"))
    }
    pheno <- pheno[missing2, , drop = FALSE]
    kinship <- kinship[missing2, missing2, drop = FALSE]
    probs <- probs[missing2, , , drop = FALSE]
    if (!is.null(addcovar)) {
        # remove subjects with missing data
        addcovar <- addcovar[missing2, , drop = FALSE]
        if (sum(!missing2) > 0){
            message(paste0("removed ",
                       sum(!missing2),
                       " subjects due to missing covariate values")
                )
        }
        # check for any addcovar that have the same value for all subjects note that we do this AFTER
        # removing subjects with missing values
        covs_identical <- apply(FUN = check_identical, X = addcovar, MARGIN = 2)
        if (sum(covs_identical) > 0) {
            message(
                paste0("removed addcovar due to absence of variation in covariate values: ",
                       colnames(addcovar)[covs_identical]))
        }
        addcovar <- addcovar[, !covs_identical, drop = FALSE]
        if (ncol(addcovar) == 0) {
            addcovar <- NULL
        }
        # remove those covariate columns for which all subjects have the same value
    }
    # create id2keep and subset all four input objects to include only those subjects that are in id2keep
    # we've already removed subjects that have missing values from both phenotypes matrix and addcovar
    # matrix
    id2keep <- make_id2keep(probs = probs, pheno = pheno, covar = addcovar, kinship = kinship)
    probs <- subset_input(input = probs, id2keep = id2keep)
    pheno <- subset_input(input = pheno, id2keep = id2keep)
    addcovar <- subset_input(input = addcovar, id2keep = id2keep)
    kinship <- subset_kinship(kinship = kinship, id2keep = id2keep)

    # covariance matrix estimation
    message("starting covariance matrices estimation.")
    # first, run gemma2::MphEM() to get Vg and Ve
    cc_out <- calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec, addcovar = addcovar)
    Vg <- cc_out$Vg
    Ve <- cc_out$Ve
    message("covariance matrices estimation completed.")

    # define Sigma
    Sigma <- calc_Sigma(Vg, Ve, kinship)
    # define Sigma_inv
    Sigma_inv <- solve(Sigma)
    # prepare table of marker indices for each call of scan_pvl
    mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
    # start progress bar
    pb <- progress::progress_bar$new(format = " scanning [:bar] :percent eta: :eta",
                                     total = n_snp^d_size,
                                     clear = FALSE,
                                     width = 80
                                     )
    pb$tick(0)

    for (rownum in 1:nrow(mytab)) {
        pb$tick()
        indices <- unlist(mytab[rownum, ])
        X_list <- prep_X_list(indices = indices[-length(indices)],
                              start_snp = start_snp,
                              probs = probs,
                              addcovar = addcovar
                              )
        X <- gemma2::stagger_mats(X_list)
        # Bhat <- rcpp_calc_Bhat2(X = X, Y = as.vector(pheno), Sigma_inv = Sigma_inv)
        Bhat <- rcpp_calc_Bhat2(X = X, Sigma_inv = Sigma_inv, Y = as.vector(as.matrix(pheno)))
        mymu <- as.vector(X %*% Bhat)
        mytab$loglik[rownum] <- as.numeric(
            rcpp_log_dmvnorm2(inv_S = Sigma_inv,
                              mu = mymu,
                              x = as.vector(as.matrix(pheno)),
                              S = Sigma)
            )
    }
    marker_id <- dimnames(probs)[[3]][start_snp:(start_snp + n_snp - 1)]
    mytab2 <- tibble::as_tibble(apply(FUN = function(x) marker_id[x], X = mytab[, -ncol(mytab)], MARGIN = 2))
    mytab2$loglik <- mytab$loglik
    return(mytab2)
}
