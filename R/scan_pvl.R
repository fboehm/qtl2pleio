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
#' @param cores number of cores to use when parallelizing via parallel::mclapply. Set to 1 for no parallelization.
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
#' Broman KW, Gatti DM, Simecek P, Furlotte NA, Prins P, Sen S, Yandell BS, Churchill GA (2019)
#' R/qtl2: software for mapping quantitative trait loci with high-dimensional data and
#' multi-parent populations. GENETICS https://www.genetics.org/content/211/2/495.
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
#' start_snp = 1, n_snp = 5, cores = 1)
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
                     cores = 1
                     )
    {
    inputs <- process_inputs(probs = probs,
                             pheno = pheno,
                             addcovar = addcovar,
                             kinship = kinship,
                             max_iter = max_iter,
                             max_prec = max_prec)
    # prepare table of marker indices for each call of scan_pvl
    d_size <- ncol(inputs$pheno)
    mytab <- prep_mytab(d_size = d_size, n_snp = n_snp)
    # set up parallel analysis
    out <- scan_pvl_clean(mytab = mytab,
                          addcovar = inputs$addcovar,
                          probs = inputs$probs,
                          Sigma_inv = inputs$Sigma_inv,
                          Sigma = inputs$Sigma,
                          start_snp = start_snp,
                          pheno = inputs$pheno,
                          n_snp = n_snp,
                          cores = cores
                          )
    return(out)
}



scan_pvl_clean <- function(pheno,
                           probs,
                           addcovar,
                           Sigma_inv,
                           Sigma,
                           start_snp,
                           mytab,
                           n_snp,
                           cores = 1){
    list_result <- parallel::mclapply(X = as.data.frame(t(mytab)),
                                     FUN = fit1_pvl,
                                     addcovar = addcovar,
                                     probs = probs,
                                     inv_S = Sigma_inv,
                                     S = Sigma,
                                     start_snp = start_snp,
                                     pheno = pheno,
                                     mc.cores = cores
                                     )
    mytab$loglik <- unlist(list_result)
    marker_id <- dimnames(probs)[[3]][start_snp:(start_snp + n_snp - 1)]
    mytab2 <- tibble::as_tibble(apply(FUN = function(x) marker_id[x],
                                      X = mytab[, -ncol(mytab)],
                                      MARGIN = 2))
    mytab2$log10lik <- mytab$loglik / log(10)
    return(mytab2)
}



