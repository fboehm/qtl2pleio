#' Perform model fitting for all ordered pairs of markers in a genomic region of interest
#'
#' @param probs a matrix of genotype probabilities for a single chromosome (before eigenrotation)
#' @param pheno a matrix of phenotypes (before eigenrotation)
#' @param kinship a kinship matrix for one chromosome
#' @param start_snp1 index of where to start the first coordinate of the ordered pair
#' @param start_snp2 index of where to start the second coordinate of the ordered pair. Typically the same value as start_snp1
#' @param n_snp the number of (consecutive) snps to include in the scan
#' @export

scan_pvl <- function(probs, pheno, kinship, start_snp1, start_snp2 = start_snp1, n_snp){
  stopifnot(nrow(probs) == nrow(pheno),
            nrow(probs) == nrow(kinship),
            nrow(kinship) == ncol(kinship)
            )

  # perform scan over probs[ , , start_snp: stop_snp]
  # first, run gemma2::MphEM() to get Vg and Ve
  gemma2::eigen2(kinship) -> e_out
  e_out$vectors -> U
  e_out$values -> eval
  n_mouse <- nrow(probs)
  X1pre <- t(rep(1, n_mouse))
  X1 <- X1pre %*% U
  Y <- t(pheno) %*% U
  # run MphEM with only a design matrix that contains only the intercept term (and not any genotype info)
  foo <- gemma2::MphEM(X = X1, Y = Y, eval = eval, V_g = diag(2), V_e = diag(2))
  Vg <- foo[[length(foo)]][[2]]
  Ve <- foo[[length(foo)]][[3]]
  # define Sigma
  Sigma <- kinship %x% Vg + diag(n_mouse) %x% Ve
  loglik <- matrix(nrow = n_snp, ncol = n_snp)
  rownames(loglik) <- dimnames(probs)[[3]][start_snp1 : (start_snp1 + n_snp - 1)]
  colnames(loglik) <- dimnames(probs)[[3]][start_snp2 : (start_snp2 + n_snp - 1)]

  for (i in 1:n_snp){
    for (j in 1:n_snp){
      index1 <- start_snp1 + i - 1
      index2 <- start_snp2 + j - 1
      X1 <- as.matrix(probs[ , , index1])
      # note that we overwrite earlier X1 here
      X2 <- as.matrix(probs[ , , index2])
      X <- pleiotropy::stagger_mats(X1, X2)
      Bhat <- calc_Bhat(X = X,
                        Y = as.vector(pheno),
                        Sigma = Sigma)
      loglik[i, j] <- calc_loglik_bvlmm(X = X, Y = pheno, Bhat = Bhat, Sigma = Sigma)
    }
  }
  return(loglik)
}
