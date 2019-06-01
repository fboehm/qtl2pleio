---
title: '`qtl2pleio`: Testing pleiotropy vs. separate QTL in multiparental populations'
tags:
  - R
  - qtl2
  - genetics
  - systems genetics
  - multiparental populations
authors:
  - name: Frederick Boehm
    orcid: 0000-0002-1644-5931
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Brian Yandell
    orcid: 0000-0002-8774-9377
    affiliation: "1"
  - name: Karl W. Broman
    orcid: 0000-0002-4914-6671
    affiliation: "2"
affiliations:
  - name: Department of Statistics, University of Wisconsin-Madison
    index: 1
  - name: Department of Biostatistics and Medical Informatics, University of Wisconsin-Madison
    index: 2
date: 26 April 2019
bibliography: paper.bib
---

# Summary

Modern quantitative trait locus (QTL) studies in multiparental populations offer opportunities to identify causal genes for thousands of clinical and molecular traits. Traditional analyses examine each trait by itself. However, to fully leverage this vast number of measured traits, the systems genetics community needs statistical tools to analyze multiple traits simultaneously [@jiang1995multiple;@korol1995interval]. A test of pleiotropy vs. separate QTL is one such tool that will aid dissection of complex trait genetics and enhance understanding of genetic architecture.

@jiang1995multiple developed a pleiotropy test for two-parent crosses. For a pair of traits that map to a single genomic region, they formulated the test with the null hypothesis being pleiotropy (the two traits are affected by a single QTL) against the alternative hypothesis of two separate QTL, with each QTL affecting exactly one trait in the pair.

The test of @jiang1995multiple doesn't directly apply to multiparental populations because

1. Multiparental populations have more than two founders
1. Multiparental populations have complicated pedigrees

Additionally, the test statistic distribution, under the null hypothesis of pleiotropy, doesn't follow a distribution with tabulated quantiles, like the chi-square with 1 degree of freedom. Thus, we need to implement a method for determining p-values. 

We addressed the first two challenges by adding a fixed effect for every founder line and incorporating a multivariate polygenic random effect into the linear model, which resulted in a multivariate linear mixed effects model [@kang2008efficient;@zhou2014efficient]. We implemented a parametric bootstrap procedure to determine p-values for test statistics [@efron1979bootstrap;@tian2016dissection]. We describe details of our statistical methods elsewhere [@boehm2019testing]. 

`qtl2pleio` offers a convenient interface for those already analyzing data with [`qtl2`](https://kbroman.org/qtl2/). The primary functions in `qtl2pleio` are `scan_pvl`, to perform the multivariate multi-QTL scan, and `boot_pvl`, to obtain bootstrap samples. We also include functions for visualizing results. `qtl2pleio` features three R package vignettes that demonstrate these and other `qtl2pleio` functions. One vignette provides examples for performing bootstrap analysis with a computing cluster. For quality assurance purposes, we incorporated unit tests into `qtl2pleio` via the R package `testthat` [@wickham2011testthat]. 

`qtl2pleio` uses `C++` code for model fitting via generalized least squares. We use the R package `Rcpp` to interface with our `C++` code [@eddelbuettel2011rcpp]. We also make use of the `C++` library `Eigen` via the R package `RcppEigen` [@bates2013fast].

# Acknowledgments

This work was supported in part by National Institutes of Health grant R01GM070683 (to K.W.B.).

# References
