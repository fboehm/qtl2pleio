#' Perform bootstrap sampling and calculate test statistics for each bootstrap sample
#'
#' Create a bootstrap sample, perform multivariate QTL scan, and calculate log10 LRT statistic
#'
#' Performs a parametric bootstrap method to calibrate test statistic values in the test of
#' pleiotropy vs. separate QTL. It begins by inferring parameter values at
#' the `pleio_peak_index` index value in the object `probs`. It then uses
#' these inferred parameter values in sampling from a multivariate normal
#' distribution. For each of the `nboot` sampled phenotype vectors, a two-dimensional QTL
#' scan, starting at the marker indexed by `start_snp` within the object
#' `probs` and extending for a total of `n_snp` consecutive markers. The
#' two-dimensional scan is performed via the function `scan_pvl_clean`. For each
#' two-dimensional scan, a log10 likelihood ratio test statistic is calculated. The
#' outputted object is a vector of `nboot` log10 likelihood ratio test
#' statistics from `nboot` distinct bootstrap samples.
#'
#' @param probs founder allele probabilities three-dimensional array for one chromosome only (not a list)
#' @param pheno n by d matrix of phenotypes
#' @param addcovar n by c matrix of additive numeric covariates
#' @param kinship a kinship matrix, not a list
#' @param start_snp positive integer indicating index within probs for start of scan
#' @param n_snp number of (consecutive) markers to use in scan
#' @param pleio_peak_index positive integer index indicating genotype matrix for bootstrap sampling. Typically acquired by using `find_pleio_peak_tib`.
#' @param nboot number of bootstrap samples to acquire and scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @export
#' @importFrom stats var
#' @references Knott SA, Haley CS (2000) Multitrait
#' least squares for quantitative trait loci detection.
#' Genetics 156: 899–911.
#'
#' Walling GA, Visscher PM, Haley CS (1998) A comparison of
#' bootstrap methods to construct confidence intervals in QTL mapping.
#' Genet. Res. 71: 171–180.
#' @examples
#'
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
#' boot_pvl(probs = probs, pheno = pheno,
#'         start_snp = 1, n_snp = 5, pleio_peak_index = 3, nboot = 1)
#'
#'
#' @return numeric vector of (log) likelihood ratio test statistics from `nboot_per_job` bootstrap samples
#'
boot_pvl <- function(probs,
                     pheno,
                     addcovar = NULL,
                     kinship = NULL,
                     start_snp = 1,
                     n_snp,
                     pleio_peak_index,
                     nboot = 1,
                     max_iter = 1e+04,
                     max_prec = 1 / 1e+08
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
    f_size <- ncol(probs)
    # force things to be matrices
    if (!is.matrix(pheno)) {
        pheno <- as.matrix(pheno)
        if (!is.numeric(pheno)) stop("pheno is not numeric")
    }
    if (is.null(colnames(pheno))) # force column names
        colnames(pheno) <- paste0("pheno", seq_len(ncol(pheno)))
    if (!is.null(addcovar)) {
        if (!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
        if (!is.numeric(addcovar)) stop("addcovar is not numeric")
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
        addcovar <- subset_input(addcovar, id2keep)
    }
    # Send messages if there are two or fewer subjects
    if (length(id2keep) == 0) {stop("no individuals common to all inputs")}
    if (length(id2keep) <= 2) {
        stop(paste0("only ", length(id2keep),
                    " common individual(s): ",
                    paste(id2keep, collapse = ": ")))
    }
    # subset inputs to get all without missingness
    probs <- subset_input(input = probs, id2keep = id2keep)
    pheno <- subset_input(input = pheno, id2keep = id2keep)
    if (!is.null(kinship)) {
        kinship <- subset_kinship(kinship = kinship, id2keep = id2keep)
    }
## define X1 - a single marker's allele probabilities
    X1 <- probs[ , , pleio_peak_index]
    if (!is.null(addcovar)) {
        Xpre <- cbind(X1, addcovar)
    } else {
        Xpre <- X1
    }
    Xlist <- vector(length = d_size)
    Xlist <- lapply(Xlist, FUN = function(x){x <- Xpre; return(x)})
    X <- gemma2::stagger_mats(Xlist)
    if (!is.null(kinship)) {
        # covariance matrix estimation
        # first, run gemma2::MphEM(), by way of calc_covs(), to get Vg and Ve
        cc_out <- calc_covs(pheno, kinship, max_iter = max_iter, max_prec = max_prec, covariates = addcovar)
        Vg <- cc_out$Vg
        Ve <- cc_out$Ve
        # define Sigma
        Sigma <- calc_Sigma(Vg = Vg, Ve = Ve, kinship = kinship)
    }
    if (is.null(kinship)) {
        # get Sigma for Haley Knott regression without random effect
        Ve <- var(pheno) # get d by d covar matrix
        Sigma <- calc_Sigma(Vg = NULL, Ve = Ve, kinship = NULL, n_mouse = nrow(pheno))
    }
    Sigma_inv <- solve(Sigma)
    # calc Bhat
    Bcol <- rcpp_calc_Bhat2(X = X,
                         Sigma_inv = Sigma_inv,
                         Y = as.vector(as.matrix(pheno))
    )
    B <- matrix(data = Bcol, nrow = ncol(Xpre), ncol = d_size, byrow = FALSE)
    # Start loop to get Ysim matrices
    Ysimlist <- list()
    for (i in 1:nboot) {
        foo <- sim1(X = X, B = B, Sigma = Sigma)
        Ysim <- matrix(foo, ncol = d_size, byrow = FALSE)
        rownames(Ysim) <- rownames(pheno)
        colnames(Ysim) <- paste0("t", 1:d_size)
        Ysimlist[[i]] <- Ysim
    }
    # prepare table of marker indices for each call of scan_pvl_clean
    mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)

    scan_out <- furrr::future_map(.x = Ysimlist,
                                   .f = scan_pvl_clean,
                                   probs = probs,
                                   addcovar = addcovar,
                                   Sigma_inv = Sigma_inv,
                                   Sigma = Sigma,
                                   start_snp = start_snp,
                                   mytab = mytab,
                                   n_snp = n_snp
                                   )
    lrt <- furrr::future_map_dbl(.x = scan_out, .f = function(x){
        x %>%
            calc_profile_lods() %>%
            dplyr::select(profile_lod) %>%
            max()
    })
    return(lrt)
}
