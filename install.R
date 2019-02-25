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
"rmarkdown")
install.packages(cran_pkgs)
devtools::install_github("rqtl/qtl2@b5785d5")
devtools::install_github("fboehm/gemma2@2872396")
devtools::install_github("fboehm/qtl2pleio@df0ffc3")
