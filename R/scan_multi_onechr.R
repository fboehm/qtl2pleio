#' Perform multivariate, one-QTL model fitting for markers on one chromosome
#'
#' `scan_multi_onechr` calculates log likelihood for d-variate phenotype model fits. Inputted parameter `start_snp` indicates where in the `probs` object to start the scan.
#'
#' @param probs an array of founder allele probabilities for a single chromosome
#' @param pheno a matrix of phenotypes
#' @param kinship a kinship matrix for one chromosome
#' @param addcovar a matrix, n subjects by c additive covariates
#' @param start_snp index of where to start the scan within probs
#' @param n_snp the number of (consecutive) markers to include in the scan
#' @param max_iter maximum number of iterations for EM algorithm
#' @param max_prec stepwise precision for EM algorithm. EM stops once incremental difference in log likelihood is less than max_prec
#' @param cores number of cores for parallelization
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
#' scan_multi_onechr(probs = probs, pheno = pheno, kinship = NULL, cores = 1)
#'
#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log10 likelihood

scan_multi_onechr <- function(probs,
                     pheno,
                     kinship = NULL,
                     addcovar = NULL,
                     start_snp = 1,
                     n_snp = dim(probs)[3],
                     max_iter = 1e+04,
                     max_prec = 1 / 1e+08,
                     cores = 1
)
{
  inputs <- process_inputs(probs = probs, pheno = pheno, kinship = kinship, addcovar = addcovar, max_iter = max_iter, max_prec = max_prec)
  # prepare table of marker indices for each call of scan_pvl_clean
  d_size <- ncol(inputs$pheno)
  mytab <- prep_mytab(d_size = d_size, n_snp = n_snp, pvl = FALSE)
  # set up parallel analysis
  out <- scan_pvl_clean(mytab = mytab,
                        addcovar = inputs$addcovar,
                        probs = inputs$probs,
                        Sigma_inv = inputs$Sigma_inv,
                        Sigma = inputs$Sigma,
                        start_snp = start_snp,
                        pheno = inputs$pheno,
                        n_snp = n_snp, cores = cores
  )
  return(out)
}

#' Perform multivariate, one-QTL model fitting for markers on all chromosomes
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
#' Generalized least squares model fitting is applied to every marker on
#' every chromosome.
#' For a single marker, we fit the model:
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

#' @param probs_list an list of arrays of founder allele probabilities
#' @param pheno a matrix of phenotypes
#' @param kinship_list a list of kinship matrices, one for each chromosome
#' @param addcovar a matrix, n subjects by c additive covariates
#' @param cores number of cores for parallelization via parallel::mclapply()
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
#' scan_multi_oneqtl(probs_list = list(probs, probs), pheno = pheno, cores = 1)
#'
#' @importFrom rlang .data
#' @return a tibble with d + 1 columns. First d columns indicate the genetic data (by listing the marker ids) used in the design matrix; last is log10 likelihood

scan_multi_oneqtl <- function(probs_list,
                              pheno,
                              kinship_list = NULL,
                              addcovar = NULL,
                              cores = 1
                              ){
  if (is.null(kinship_list)) {out_list <- parallel::mclapply(X = probs_list,
                                 mc.cores = cores,
                                 FUN = function(x){
                                                  scan_multi_onechr(probs = x,
                                                                    pheno = pheno,
                                                                    kinship = NULL,
                                                                    addcovar = addcovar, cores = cores
                                                  )
                                 })} else {
                                   out_list <- parallel::mclapply(X = probs_list,
                                                                  mc.cores = cores,
                                                                 FUN = function(x){
                                                                   scan_multi_onechr(probs = x,
                                                                                     pheno = pheno,
                                                                                     kinship = kinship_list,
                                                                                     addcovar = addcovar, cores = cores
                                                                   )
                                                                 })
                                                }
  # calculate log10lik for "null" model without genotypes
  if (!is.null(kinship_list)){
    inputs <- process_inputs(probs = probs_list[[1]], # arbitrary choice of which array
                           pheno = pheno,
                           addcovar = addcovar,
                           kinship = kinship_list[[1]] #arbitrary choice of which kin matrix
                           )
  } else {
    inputs <- process_inputs(probs = probs_list[[1]], # arbitrary choice of which array
                             pheno = pheno,
                             addcovar = addcovar,
                             kinship = kinship_list #arbitrary choice of which kin matrix
    )
  }
  n <- nrow(inputs$pheno)
  d_size <- ncol(inputs$pheno)
  if (is.null(addcovar)) {
    Xlist <- lapply(X = as.list(1:d_size),
                  FUN = function(x){return(matrix(data = rep(1, n),
                                                               nrow = n,
                                                               ncol = 1))})
  } else {
    Xlist <- lapply(X = as.list(1:d_size),
                          FUN = function(x){return(cbind(inputs$addcovar, matrix(data = rep(1, n),
                                                          nrow = n,
                                                          ncol = 1)))})
  }
  X <- gemma2::stagger_mats(Xlist)
  Bhat <- rcpp_calc_Bhat2(X = X,
                          Sigma_inv = inputs$Sigma_inv,
                          Y = as.vector(as.matrix(inputs$pheno))
  )
  mymu <- as.vector(X %*% Bhat)
  out <- rcpp_log_dmvnorm2(inv_S = inputs$Sigma_inv,
                           mu = mymu,
                           x = as.vector(as.matrix(inputs$pheno)),
                           S = inputs$Sigma
                           ) / log(10)
  out_list <- parallel::mclapply(X = out_list, mc.cores = cores,
                    FUN = function(x){
                      x %>%
                        dplyr::mutate(null_log10lik = as.numeric(out))

  })

  return(out_list)
  }

#' Permute the phenotypes matrix and then scan the genome. Record the genomewide greatest LOD score for each permuted data set.
#'
#' @param probs_list a list of founder allele probabilities, one array per chromosome
#' @param pheno a matrix of trait values
#' @param kinship_list a list of kinship matrices, one per chromosome
#' @param addcovar a matrix of covariate values
#' @param n_perm positive integer for the number of permuted data sets to scan.
#' @param cores number of cores for parallelization
#' @return a vector of `n_perm` max lod statistics
#' @export

scan_multi_oneqtl_perm <- function(probs_list,
                              pheno,
                              kinship_list = NULL,
                              addcovar = NULL,
                              n_perm = 1,
                              cores = 1
){
  # create a permuted phenotypes matrix
  n <- nrow(pheno)
  phe <- list() # list of permuted data sets
  for (i in 1:n_perm){
    permuted_pheno <- pheno[sample(1:n, replace = FALSE), , drop = FALSE]
    rownames(permuted_pheno) <- rownames(pheno)
    phe[[i]] <- permuted_pheno
  }
  parallel::mclapply(X = phe, mc.cores = cores, FUN = function(x){
    sout <- scan_multi_oneqtl(probs_list = probs_list,
                      pheno = x,
                      kinship_list = kinship_list,
                      addcovar = addcovar
                      )
    maxlod_perchr <- parallel::mclapply(X = sout, mc.cores = cores, FUN = function(x){
      x %>%
        dplyr::select(lod) %>%
        max()
    })
    maxlod <- max(maxlod_perchr)
    return(maxlod)
  })
}
