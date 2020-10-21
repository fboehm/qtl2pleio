#' Process inputs to scan functions
#'
#' @param probs a three-dimensional array of founder allele probabilities
#' @param pheno a matrix of d trait values
#' @param addcovar a matrix of covariates
#' @param kinship a kinship matrix
#' @param n_snp number of markers
#' @param start_snp index number of start position in the probs object.
#' @param max_iter max number of iterations for EM
#' @param max_prec max precision for stopping EM

process_inputs <- function(probs,
                           pheno,
                           addcovar,
                           kinship,
                           n_snp = dim(probs)[3],
                           start_snp = 1,
                           max_iter = 10 ^ 4,
                           max_prec = 1 / 10 ^ 8){
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
    addcovar <- subset_input(input = addcovar, id2keep = id2keep)
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
  if (d_size != qr(pheno)$rank) stop("Phenotypes matrix is not full rank. Input only full-rank phenotypes matrices.")

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
  # make list to return
  out <- list(probs = probs,
              pheno = pheno,
              addcovar = addcovar,
              kinship = kinship,
              Sigma = Sigma,
              Sigma_inv = Sigma_inv)
  return(out)
}
