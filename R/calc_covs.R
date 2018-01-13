#' Calculate Vg and Ve from Y and kinship
#'
#' @param pheno n by d matrix of phenotypes
#' @param kinship a kinship matrix, n by n
#' @param X1pre n by c design matrix. c = 1 to ignore genotypes
#' @export
calc_covs <- function(pheno, kinship, X1pre = rep(1, nrow(kinship))){
  gemma2::eigen2(kinship) -> e_out
  e_out$vectors -> U
  e_out$values -> eval
  n_mouse <- nrow(kinship)
  X1 <- t(X1pre) %*% U
  Y <- t(pheno) %*% U
  ncol(pheno) -> d
  # run MphEM with only a design matrix that contains only the intercept term (and not any genotype info)
  foo <- gemma2::MphEM(X = X1, Y = Y, eval = eval, V_g = diag(d), V_e = diag(d))
  Vg <- foo[[length(foo)]][[2]]
  Ve <- foo[[length(foo)]][[3]]
  out <- list(Vg = Vg, Ve = Ve)
  return(out)
}
