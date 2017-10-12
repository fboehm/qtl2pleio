#' Calculate Vg and Ve from Y and kinship
#'
#' @param pheno n by 2 matrix of phenotypes
#' @param kinship a kinship matrix, n by n
#' @export
calc_covs <- function(pheno, kinship){
  gemma2::eigen2(kinship) -> e_out
  e_out$vectors -> U
  e_out$values -> eval
  n_mouse <- nrow(kinship)
  X1pre <- t(rep(1, n_mouse))
  X1 <- X1pre %*% U
  Y <- t(pheno) %*% U
  # run MphEM with only a design matrix that contains only the intercept term (and not any genotype info)
  foo <- gemma2::MphEM(X = X1, Y = Y, eval = eval, V_g = diag(2), V_e = diag(2))
  Vg <- foo[[2]]
  Ve <- foo[[3]]
  out <- list(Vg = Vg, Ve = Ve)
  return(out)
}
