#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' @param probs an array of genotype probabilities for a single chromosome (before eigenrotation)
#' @param pheno a matrix of phenotypes (before eigenrotation)
#' @param kinship a kinship matrix for one chromosome
#' @param covariates a matrix, n subjects by n.cov covariates, where each column is one covariate
#' @param start_snp1 index of where to start the first coordinate of the ordered pair
#' @param start_snp2 index of where to start the second coordinate of the ordered pair. Typically the same value as start_snp1
#' @param n_snp the number of (consecutive) snps to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @export

scan_pvl <- function(probs, pheno, kinship, covariates = NULL, start_snp1,
                     start_snp2 = start_snp1, n_snp, max_iter = 100000,
                     max_prec = 1 / 1e06){
  stopifnot(identical(nrow(probs), nrow(pheno)), identical(rownames(probs), rownames(pheno)),
            identical(rownames(kinship), rownames(pheno)),
            n_snp > 0,
            start_snp1 > 0,
            start_snp1 + n_snp - 1 <= dim(probs)[3]
            )
  # start progress bar
  pb <- progress::progress_bar$new(
    format = " scanning [:bar] :percent eta: :eta",
    total = n_snp * n_snp, clear = FALSE, width= 80)
  pb$tick(0)
  # remove mice with missing values of phenotype or missing value(s) in covariates
  missing_indic <- matrix(!apply(FUN = is.finite, X = pheno, MARGIN = 1),
                          nrow = nrow(pheno), ncol = ncol(pheno),
                          byrow = TRUE)
  missing2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(pheno))), MARGIN = 1, X = missing_indic)
  if (!is.null(covariates)){
    miss_cov <- matrix(!apply(FUN = is.finite, X = covariates, MARGIN = 1),
                       nrow = nrow(covariates), ncol = ncol(covariates),
                       byrow = TRUE)
    miss_cov2 <- apply(FUN = function(x)identical(as.logical(x), rep(FALSE, ncol(covariates))), MARGIN = 1, X = miss_cov)
    missing2 <- missing2 & miss_cov2
  }
  if (sum(!missing2) > 0){message(paste0(sum(!missing2), " subjects dropped due to missing values"))}
  pheno2 <- pheno[missing2, , drop = FALSE]
  pheno2 -> pheno
  kinship <- kinship[missing2, missing2, drop = FALSE]
  probs <- probs[missing2, , , drop = FALSE]
  # perform scan over probs[ , , start_snp: stop_snp]
  # first, run gemma2::MphEM() to get Vg and Ve
  calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec) -> cc_out
  Vg <- cc_out$Vg
  Ve <- cc_out$Ve
  # define Sigma
  n_mouse <- nrow(kinship)
  Sigma <- calc_Sigma(Vg, Ve, kinship)
  # define Sigma_inv
  Sigma_inv <- solve(Sigma)
  loglik <- matrix(nrow = n_snp, ncol = n_snp)
  rownames(loglik) <- dimnames(probs)[[3]][start_snp1 : (start_snp1 + n_snp - 1)]
  colnames(loglik) <- dimnames(probs)[[3]][start_snp2 : (start_snp2 + n_snp - 1)]

  for (i in 1:n_snp){
    for (j in 1:n_snp){
      pb$tick()
      index1 <- start_snp1 + i - 1
      index2 <- start_snp2 + j - 1
      if (!is.null(covariates)){
        X1 <- cbind(as.matrix(probs[ , , index1]), covariates[missing2, , drop = FALSE])
        # note that we overwrite earlier X1 here
        X2 <- cbind(as.matrix(probs[ , , index2]), covariates[missing2, , drop = FALSE])
      } else {
        X1 <- as.matrix(probs[ , , index1])
        # note that we overwrite earlier X1 here
        X2 <- as.matrix(probs[ , , index2])
      }
      X <- gemma2::stagger_mats(X1, X2)
      Bhat <- calc_Bhat(X = X,
                        Y = as.vector(pheno),
                        Sigma_inv = Sigma_inv)
      loglik[i, j] <- calc_loglik_bvlmm(X = X, Y = as.vector(as.matrix(pheno)), Bhat = Bhat, Sigma = Sigma)
    }
  }
  return(loglik)
}
