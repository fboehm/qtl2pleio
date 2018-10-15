#' Create a bootstrap sample, perform multivariate QTL scan, and calculate LRT statistic
#'
#'
#'
#'
#'
#' @param probs allele probabilities object for one chromosome only, like aprobs$`8`. Not a list
#' @param pleio_peak_index positive integer index indicating design matrix for simulation. Typically acquired by using `find_pleio_peak_tib`.
#' @param pheno n by d matrix of phenotypes
#' @param covariates n by n.cov matrix of covariates
#' @param kinship a kinship matrix, not a list
#' @param nboot_per_job number of bootstrap samples to call per invocation of function
#' @param start_snp positive integer indicating index within probs for start of scan
#' @param n_snp number of (consecutive) markers to use in scan
#' @export
#' @return numerical vector of lrt statistics from `nboot_per_job` bootstrap samples
#'
boot_pvl <- function(probs,
                     pleio_peak_index,
                     pheno,
                     covariates = NULL,
                     kinship = NULL,
                     nboot_per_job = 1,
                     start_snp,
                     n_snp) {
    stopifnot(identical(nrow(kinship), nrow(probs)),
              identical(nrow(probs), nrow(pheno)),
              check_dimnames(kinship, probs),
              check_dimnames(probs, pheno))
    if (!is.null(covariates)) {
        stopifnot(check_dimnames(kinship, covariates),
                  identical(nrow(probs), nrow(covariates))
                  )
    }
    X1 <- probs[, , pleio_peak_index]
    if (!is.null(covariates)) {
        Xpre <- cbind(X1, covariates)
    } else {
        Xpre <- X1
    }
    ## remove subjects with missing values of phenotype
    missing_indic <- apply(FUN = function(x) identical(x, rep(TRUE, length(x))), X = !is.na(pheno),
        MARGIN = 1)
    pheno_nona <- pheno[!missing_indic, ]
    Xpre_nona <- Xpre[!missing_indic, ]
    k_nona <- kinship[!missing_indic, !missing_indic]
    ##
    X <- gemma2::stagger_mats(Xpre_nona, Xpre_nona)
    cc_out <- calc_covs(pheno = pheno_nona,
                        kinship = k_nona,
                        covariates = covariates
                        )
    (Vg <- cc_out$Vg)
    (Ve <- cc_out$Ve)
    # calculate Sigma
    Sigma <- calc_Sigma(Vg = Vg,
                        Ve = Ve,
                        K = k_nona
                        )
    Sigma_inv <- solve(Sigma)
    # calc Bhat
    B <- calc_Bhat(X = X,
                   Sigma_inv = Sigma_inv,
                   Y = pheno_nona
                   )
    # Start loop
    lrt <- numeric()
    for (i in 1:nboot_per_job) {
        foo <- sim1(X = X, B = B, Vg = Vg, Ve = Ve, kinship = k_nona)
        Ysim <- matrix(foo, ncol = 2, byrow = FALSE)
        rownames(Ysim) <- rownames(pheno_nona)
        colnames(Ysim) <- c("t1", "t2")
        if (!is.null(covariates)) {
            loglik <- scan_pvl(probs = probs[!missing_indic, , ],
                               pheno = Ysim,
                               covariates = covariates[!missing_indic, , drop = FALSE],
                               kinship = k_nona,
                               start_snp = start_snp,
                               n_snp = n_snp
                               )
        } else {
            loglik <- scan_pvl(probs = probs[!missing_indic, , ],
                               pheno = Ysim,
                               kinship = k_nona,
                               start_snp = start_snp,
                               n_snp = n_snp
                               )
        }
        lrt[i] <- calc_lrt_tib(loglik)
    }
    return(lrt)
}
