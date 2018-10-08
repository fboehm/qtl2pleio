#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' `scan_pvl` calculates log likelihood for d-variate phenotype model fits. Inputted parameter `start_snp1` indicates where in the `probs` object to start the scan.
#'
#' The function first discards individuals with missing phenotypes or missing covariates.
#' It also discards covariates that have the same value for all remaining subjects.
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
#' @param covariates a matrix, n subjects by n.cov covariates, where each column is one covariate
#' @param start_snp1 index of where to start the scan within probs
#' @param n_snp the number of (consecutive) markers to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @param use_limmbo2 logical indicating whether to use limmbo2::limbbo2() in covariance matrices estimation. Default is to use gemma2 R package.
#' @param limmbo2_subset_size positive integer denoting the number of phenotypes for each bootstrap sample with limmbo2. Only used when `use_limmbo2` is TRUE
#' @export
#' @examples
#' ## define probs
#'probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
#'probs <- array(data = probs_pre, dim = c(100, 1, 10))
#'s_id <- paste0("s", 1:100)
#'rownames(probs) <- s_id
#'colnames(probs) <- "A"
#'dimnames(probs)[[3]] <- paste0("Marker", 1:10)
#'# define Y
#'Y_pre <- runif(200)
#'Y <- matrix(data = Y_pre, nrow = 100)
#'rownames(Y) <- s_id
#'colnames(Y) <- paste0("t", 1:2)
#'covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
#'rownames(covariates) <- s_id
#'colnames(covariates) <- "c1"
#'Y2 <- Y
#'Y2[1, 2] <- NA
#'scan_pvl(probs = probs, pheno = Y, kinship = diag(100),
#'         start_snp1 = 1, n_snp = 10)
#'scan_pvl(probs = probs, pheno = Y2, kinship = diag(100),
#'         start_snp1 = 1, n_snp = 10)
#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log likelihood

scan_pvl <- function(probs, pheno, kinship,
                     covariates = NULL,
                     start_snp1 = 1,
                     n_snp,
                     max_iter = 1000000,
                     max_prec = 1 / 1e08,
                     use_limmbo2 = FALSE,
                     limmbo2_subset_size = NULL){
  stopifnot(
            !is.null(rownames(probs)),
            !is.null(colnames(probs)),
            !is.null(dimnames(probs)[[3]]),
            !is.null(rownames(pheno)),
            !is.null(colnames(pheno)),
            !is.null(rownames(kinship)),
            !is.null(colnames(kinship)),
            n_snp > 0,
            start_snp1 > 0,
            start_snp1 + n_snp - 1 <= dim(probs)[3]
            )
  # check additional conditions when covariates is not NULL
  if(!is.null(covariates)){
    stopifnot(
      !is.null(rownames(covariates)),
      !is.null(colnames(covariates))
      )}

  d_size <- ncol(pheno) # d_size is the number of univariate phenotypes
  # remove mice with missing values of phenotype or missing value(s) in covariates
  missing_indic <- matrix(!apply(FUN = is.finite, X = pheno, MARGIN = 1),
                          nrow = nrow(pheno), ncol = ncol(pheno),
                          byrow = TRUE)
  missing2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(pheno))),
                    MARGIN = 1, X = missing_indic)
  if (!is.null(covariates)){
    covariates <- subset_input(input = covariates, id2keep = intersect(rownames(pheno), rownames(covariates)))
    pheno <- subset_input(input = pheno, id2keep = intersect(rownames(pheno), rownames(covariates)))
    miss_cov <- matrix(!apply(FUN = is.finite, X = covariates, MARGIN = 1),
                       nrow = nrow(covariates), ncol = ncol(covariates),
                       byrow = TRUE)
    miss_cov2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(covariates))),
                       MARGIN = 1, X = miss_cov)
    missing2 <- missing2 & miss_cov2
  }
  if (sum(!missing2) > 0){
    message(paste0(sum(!missing2), " subjects dropped due to missing values"))
  }
  pheno <- pheno[missing2, , drop = FALSE]
  kinship <- kinship[missing2, missing2, drop = FALSE]
  probs <- probs[missing2, , , drop = FALSE]
  if (!is.null(covariates)){
    # remove subjects with missing data
    covariates <- covariates[missing2, , drop = FALSE]
    # check for any covariates that have the same value for all subjects
    # note that we do this AFTER removing subjects with missing values
    covs_identical <- apply(FUN = check_identical, X = covariates, MARGIN = 2)
    covariates <- covariates[ , !covs_identical, drop = FALSE]
    if (ncol(covariates) == 0){
      covariates <- NULL
    }
    if (sum(covs_identical) > 0){
      message(paste0(sum(covs_identical), " covariates dropped due to no variation in covariates"))
    }
    # remove those covariate columns for which all subjects have the same value
  }
  # create id2keep and subset all four input objects to include only those subjects that are in id2keep
  # we've already removed subjects that have missing values from both phenotypes matrix and covariates matrix
  id2keep <- make_id2keep(probs = probs,
                          pheno = pheno,
                          covar = covariates,
                          kinship = kinship)
  probs <- subset_input(input = probs, id2keep = id2keep)
  pheno <- subset_input(input = pheno, id2keep = id2keep)
  covariates <- subset_input(input = covariates, id2keep = id2keep)
  kinship <- subset_kinship(kinship = kinship, id2keep = id2keep)

  # covariance matrix estimation
  message("starting covariance matrices estimation.")
  if (!use_limmbo2){
    # first, run gemma2::MphEM() to get Vg and Ve
    calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec,
              covariates = covariates) -> cc_out
    Vg <- cc_out$Vg
    Ve <- cc_out$Ve
  }
  else {
    limmbo2::limmbo2(kinship = kinship, pheno = pheno, S = limmbo2_subset_size) -> li_out
    Vg <- li_out$Vg
    Ve <- li_out$Ve
  }
  message("covariance matrices estimation completed.")

  # define Sigma
  n_mouse <- nrow(kinship) # define n_mouse as the number of mice that actually have no missing data
  Sigma <- calc_Sigma(Vg, Ve, kinship)
  # define Sigma_inv
  Sigma_inv <- solve(Sigma)
  mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
  # start progress bar
  pb <- progress::progress_bar$new(
    format = " scanning [:bar] :percent eta: :eta",
    total = n_snp ^ d_size, clear = FALSE, width = 80)
  pb$tick(0)

  for (rownum in 1:nrow(mytab)){
    pb$tick()
    indices <- unlist(mytab[rownum, ])
    X_list <- prep_X_list(indices = indices[ - length(indices)],
                          start_snp1 = start_snp1,
                          probs = probs,
                          covariates = covariates)
    X <- gemma2::stagger_mats(X_list)
    #Bhat <- rcpp_calc_Bhat2(X = X,
    #                  Y = as.vector(pheno),
    #                  Sigma_inv = Sigma_inv)
    Bhat <- rcpp_calc_Bhat(X = X, Sigma = Sigma, Y = as.vector(as.matrix(pheno)))
    as.vector(X %*% Bhat) -> mymu
    mytab$loglik[rownum] <- as.numeric(rcpp_log_dmvnorm2(inv_S = Sigma_inv, mu = mymu,
                                             x = as.vector(as.matrix(pheno)),
                                             S = Sigma
                                             ))
  }
  marker_id <- dimnames(probs)[[3]][start_snp1:(start_snp1 + n_snp - 1)]
  tibble::as_tibble(apply(FUN = function(x)marker_id[x],
                          X = mytab[ , - ncol(mytab)], MARGIN = 2)) -> mytab2
  mytab2$loglik <- mytab$loglik
  return(mytab2)
}
