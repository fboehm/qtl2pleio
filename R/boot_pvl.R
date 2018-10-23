#' Perform bootstrap sampling and calculate test statistics for each bootstrap sample
#'
#' Create a bootstrap sample, perform multivariate QTL scan, and calculate LRT statistic
#'
#' Performs a parametric bootstrap method to calibrate test statistic values in the test of
#' pleiotropy vs. separate QTL. It begins by inferring parameter values at
#' the `pleio_peak_index` index value in the object `probs`. It then uses
#' these inferred parameter values in sampling from a multivariate normal
#' distribution. For each of the `nboot_per_job` sampled phenotype vectors, a two-dimensional QTL
#' scan, starting at the marker indexed by `start_snp` within the object
#' `probs` and extending for a total of `n_snp` consecutive markers. The
#' two-dimensional scan is performed via the function `scan_pvl`. For each
#' two-dimensional scan, a likelihood ratio test statistic is calculated. The
#' outputted object is a vector of `nboot_per_job` likelihood ratio test
#' statistics from `nboot_per_job` distinct bootstrap samples.
#'
#'
#'
#'
#' @param probs allele probabilities object for one chromosome only, like aprobs$`8`. Not a list
#' @param pheno n by d matrix of phenotypes
#' @param addcovar n by n.cov matrix of additive numeric covariates
#' @param kinship a kinship matrix, not a list
#' @param start_snp positive integer indicating index within probs for start of scan
#' @param n_snp number of (consecutive) markers to use in scan
#' @param pleio_peak_index positive integer index indicating design matrix for simulation. Typically acquired by using `find_pleio_peak_tib`.
#' @param nboot_per_job number of bootstrap samples to call per invocation of function
#' @export
#' @references Knott SA, Haley CS (2000) Multitrait
#' least squares for quantitative trait loci detection.
#' Genetics 156: 899–911.
#'
#' Walling GA, Visscher PM, Haley CS (1998) A comparison of
#' bootstrap methods to construct confidence intervals in QTL mapping.
#' Genet. Res. 71: 171–180.
#' @examples
#'
## define probs
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
#'addcovar <- matrix(c(runif(99), NA), nrow = 100, ncol = 1)
#'rownames(addcovar) <- s_id
#'colnames(addcovar) <- 'c1'
#'kin <- diag(100)
#'rownames(kin) <- s_id
#'colnames(kin) <- s_id
#'Y2 <- Y
#'Y2[1, 2] <- NA
#'set.seed(2018-10-22)
#'boot_pvl(probs = probs, pheno = Y, kinship = kin,
#'         start_snp = 1, n_snp = 10, pleio_peak_index = 10, nboot_per_job = 1)
#'boot_pvl(probs = probs, pheno = Y2, kinship = kin,
#'         start_snp = 1, n_snp = 10, pleio_peak_index = 10, nboot_per_job = 2)
#'
#'
#' @return numeric vector of lrt statistics from `nboot_per_job` bootstrap samples
#'
boot_pvl <- function(probs,
                     pheno,
                     addcovar = NULL,
                     kinship = NULL,
                     start_snp = 1,
                     n_snp,
                     pleio_peak_index,
                     nboot_per_job = 1
                     )
    {
    if(is.null(probs)) stop("probs is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    stopifnot(identical(nrow(kinship), nrow(probs)),
              identical(nrow(probs), nrow(pheno)),
              check_dimnames(kinship, probs),
              check_dimnames(probs, pheno))
    if (!is.null(addcovar)) {
        stopifnot(check_dimnames(kinship, addcovar),
                  identical(nrow(probs), nrow(addcovar))
                  )
    }




    X1 <- probs[, , pleio_peak_index]
    if (!is.null(addcovar)) {
        Xpre <- cbind(X1, addcovar)
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
                        covariates = addcovar
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
        if (!is.null(addcovar)) {
            loglik <- scan_pvl(probs = probs[!missing_indic, , ],
                               pheno = Ysim,
                               addcovar = addcovar[!missing_indic, , drop = FALSE],
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
