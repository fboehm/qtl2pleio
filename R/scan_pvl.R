#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' @param probs an array of genotype probabilities for a single chromosome (before eigenrotation)
#' @param pheno a matrix of phenotypes (before eigenrotation)
#' @param kinship a kinship matrix for one chromosome
#' @param covariates a matrix, n subjects by n.cov covariates, where each column is one covariate
#' @param start_snp1 index of where to start the scan within probs
#' @param n_snp the number of (consecutive) snps to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @param use_limmbo2 logical indicating whether to use limmbo2::limbbo2() in covariance matrices estimation. Default is to use gemma2 R package.
#' @param limmbo2_subset_size positive integer denoting the number of phenotypes for each bootstrap sample with limmbo2. Only used when `use_limmbo2` is TRUE
#' @export
#' @examples
#' ## define probs
#' probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
#'probs <- array(data = probs_pre, dim = c(100, 1, 10))
#'dimnames(probs)[[3]] <- paste0("Marker", 1:10)
#'# define Y
#'Y_pre <- runif(200)
#'Y <- matrix(data = Y_pre, nrow = 100)
#'covariates <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
#'cov2 <- matrix(c(covariates[1:99], 10), nrow = 100, ncol = 1)
#'Y2 <- Y
#'Y2[1, 2] <- NA
#'scan_pvl(probs = probs, pheno = Y, kinship = diag(100), start_snp1 = 1, n_snp = 10)
#'scan_pvl(probs = probs, pheno = Y2, kinship = diag(100), start_snp1 = 1, n_snp = 10)

#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log likelihood

scan_pvl <- function(probs, pheno, kinship, covariates = NULL, start_snp1 = 1,
                     n_snp, max_iter = 1000000,
                     max_prec = 1 / 1e08, use_limmbo2 = FALSE,
                     limmbo2_subset_size = NULL){
  stopifnot(identical(nrow(probs), nrow(pheno)),
            identical(rownames(probs), rownames(pheno)),
            identical(rownames(kinship), rownames(pheno)),
            n_snp > 0,
            start_snp1 > 0,
            start_snp1 + n_snp - 1 <= dim(probs)[3],
            if(!is.null(covariates)){identical(rownames(covariates), rownames(probs))}
            )
  d_size <- ncol(pheno)
  # start progress bar
  pb <- progress::progress_bar$new(
    format = " scanning [:bar] :percent eta: :eta",
    total = n_snp ^ d_size, clear = FALSE, width= 80)
  pb$tick(0)
  ## define number of dimensions, d_size
  # remove mice with missing values of phenotype or missing value(s) in covariates
  missing_indic <- matrix(!apply(FUN = is.finite, X = pheno, MARGIN = 1),
                          nrow = nrow(pheno), ncol = ncol(pheno),
                          byrow = TRUE)
  missing2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(pheno))),
                    MARGIN = 1, X = missing_indic)
  if (!is.null(covariates)){
    miss_cov <- matrix(!apply(FUN = is.finite, X = covariates, MARGIN = 1),
                       nrow = nrow(covariates), ncol = ncol(covariates),
                       byrow = TRUE)
    miss_cov2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(covariates))),
                       MARGIN = 1, X = miss_cov)
    missing2 <- missing2 & miss_cov2
  }
  if (sum(!missing2) > 0){message(paste0(sum(!missing2), " subjects dropped due to missing values"))}
  pheno <- pheno[missing2, , drop = FALSE]
  kinship <- kinship[missing2, missing2, drop = FALSE]
  probs <- probs[missing2, , , drop = FALSE]
  if (!is.null(covariates)){
    # remove subjects with missing data
    covariates <- covariates[missing2, , drop = FALSE]
    # check for any covariates that have the same value for all subjects
    # note that we do this AFTER removing subjects with missing values
    covs_identical <- apply(FUN = check_identical, X = covariates, MARGIN = 2)
    covariates <- covariates[ , !covs_identical]
    if(ncol(covariates) == 0){covariates <- NULL}
    if (sum(covs_identical) > 0){message(paste0(sum(covs_identical), " covariates dropped due to no variation in covariates"))}
    # remove those covariate columns
    # for which all subjects have the same value
  }

  # perform scan over probs[ , , start_snp: stop_snp]
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
  # define Sigma
  n_mouse <- nrow(kinship)
  Sigma <- calc_Sigma(Vg, Ve, kinship)
  # define Sigma_inv
  Sigma_inv <- solve(Sigma)
  mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
  for (rownum in 1:nrow(mytab)){
    pb$tick()
    indices <- unlist(mytab[rownum, ])
    X_list <- prep_X_list(indices = indices[ - length(indices)],
                          start_snp1 = start_snp1, probs = probs,
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
