% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_Bhat.R
\name{calc_Bhat}
\alias{calc_Bhat}
\title{Calculate estimated allele effects, B matrix}
\usage{
calc_Bhat(X, Sigma_inv, Y)
}
\arguments{
\item{X}{dn by df block-diagonal design matrix that incorporates genetic info for d markers. Note that we can use the same marker data twice.}

\item{Sigma_inv}{dn by dn inverse covariance matrix, often composed as the inverse of \eqn{K \otimes V_g + I_n \otimes V_e}}

\item{Y}{dn by 1 matrix, ie, a column vector, of d phenotypes' measurements}
}
\value{
a df by 1 matrix of GLS-estimated allele effects
}
\description{
Calculate estimated allele effects, B matrix
}
\examples{
X1 <- as.matrix(rbinom(n = 100, size = 1, prob = 1 / 2))
X <- gemma2::stagger_mats(X1, X1)
Sigma_inv <- diag(200)
Y <- runif(200)
calc_Bhat(X, Sigma_inv, Y)
}
