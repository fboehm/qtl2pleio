cran_pkgs<- c("MASS",
"microbenchmark",
"mvtnorm",
"progress",
"Rcpp",
"RcppEigen",
"rlang",
"parallel",
"broman",
"devtools",
"rmarkdown", "qtl2", "gemma2", "furrr")
install.packages("versions")
versions::install.dates(cran_pkgs, dates = "2020-07-06")
