% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_identical}
\alias{check_identical}
\title{Check whether a vector, x, has all its entries equal to its first entry}
\usage{
check_identical(x)
}
\arguments{
\item{x}{a vector}
}
\value{
a logical indicating whether all vector entries are the same
}
\description{
Check whether a vector, x, has all its entries equal to its first entry
}
\examples{
x <- 1:5
check_identical(x)
y <- rep(1, 5)
check_identical(y)
}
