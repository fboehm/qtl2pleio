% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_mytab.R
\name{prep_mytab}
\alias{prep_mytab}
\title{Prepare mytab object for use within scan_pvl R code}
\usage{
prep_mytab(d_size, n_snp, pvl = TRUE)
}
\arguments{
\item{d_size}{an integer, the number of traits}

\item{n_snp}{an integer, the number of markers}

\item{pvl}{logical indicating whether to output dataframe with all d-tuples for a d-QTL scan, or only those models that examine one marker at a time.}
}
\value{
a data.frame with d_size + 1 columns and (n_snp)^d_size rows. Last column is NA and named loglik.
}
\description{
Prepare mytab object for use within scan_pvl R code
}
\examples{
prep_mytab(2, 10)
}
