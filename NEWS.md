# qtl2pleio 1.4.1

## Major changes
 
* None

## Minor changes

* Added vignette to .Rbuildignore to ease CRAN submission.


# qtl2pleio 1.4.0

## Major changes

* Refactored `scan_pvl` and `boot_pvl` to reduce duplicated code
* Added functions for genomewide, multivariate, one-QTL scan

## Minor changes

* Returned to parallelization implemented by R package `parallel`.
* Set default `cores` value to `parallelly::availableCores()` for compatibility with job schedulers.



# qtl2pleio 1.3.0

## Major changes

* Replaced `parallel` package use with `furrr`. This led to removal of the argument 
`n_cores` from both `boot_pvl` and `scan_pvl`. A user may now use `future::plan` to 
run in parallel.
* Moved vignettes to website: https://fboehm.us/software/qtl2pleio

## Minor changes

* Fixed argument name in `boot_pvl`, changing `nboot_per_job` to `nboot`.
* Created function `scan_pvl_clean` to be called in both `scan_pvl` and `boot_pvl`. This change shouldn't affect most users' experiences.

# qtl2pleio 1.2.4

## Major changes 

* None

## Minor changes

* Updated `calc_lrt_tib` function to accommodate d-dimensional, d-variate scans


# qtl2pleio 1.2.3

## Major changes

* None

## Minor changes

* I shortened title per CRAN recommendation

* I added \value tags to document functions' outputs


# qtl2pleio 1.2.2

## Major changes

* I added data to package to enable vignettes, examples, and tests.

## Minor changes

* Per CRAN submission feedback, I spelled out QTL in DESCRIPTION file and upon first use in vignettes.

* Added citations with references to vignettes.


# qtl2pleio 1.2.1

## Major changes

* none

## Minor changes

* removed examples' and vignettes' dependencies on qtl2 R package.


# qtl2pleio 1.2.0
 
## Major changes

* changed output of scan_pvl to be log10 likelihood, not log (natural base) likelihood. Thank you to @HongHe0123 for suggesting this.

## Minor changes

* changed column name for tibble outputs to reflect log base 10 instead of log natural base. 



# qtl2pleio 1.1.0

## Major changes 

* rewrites (and deletions) of tibble helper functions
* major revisions to code in main vignette

## Minor changes

* fixed CITATION after G3 publication

# qtl2pleio 1.0.7

## Major changes

* None

## Minor changes

* fixed Description field of DESCRIPTION file per feedback from CRAN

# qtl2pleio 1.0.6

## Major changes

* None

## Minor changes

* removed Additional_repositories field from DESCRIPTION


# qtl2pleio 1.0.5

## Major changes

* None

## Minor changes

* added cran-comments.md file via devtools & added it to .Rbuildignore

# qtl2pleio 1.0.4

## Major changes

* None

## Minor changes

* Added Additional_repositories field to DESCRIPTION file per https://cran.r-project.org/submit.html in preparation for CRAN submission.  
* Added url for Jiang and Zeng 1995 GENETICS article, per https://cran.r-project.org/submit.html.

# qtl2pleio 1.0.3

## Major changes

* Revisions per peer review at Journal of Open Source Software

## Minor changes

* None

# qtl2pleio 0.1.3

## Major changes

* Submitted paper.md to JOSS

## Minor changes

* Minor updates to recla vignette to download precalculated founder allele dosages


# qtl2pleio 0.1.2.9002

## Major changes

* changed default behavior of `plot_pvl` to not indicate univariate peak positions. 
* fixed CITATION file to tell users to cite Biorxiv preprint (following CITATION file format of rqtl/qtl2)

## Minor changes

* Removed hex logo file from repo and deleted it from README
* revised paper.md in preparation for JOSS submission


## Bug fixes

* Re-created profile LOD plots for vignettes and README without indicating univariate peak positions.

* Added line type arguments to `plot_pvl`.


# qtl2pleio 0.1.2.9001

## Major changes

* added a vignette for HTCondor & bootstrap analysis
* updated README.Rmd and README.md (per suggestions & PR of @kbroman)


# qtl2pleio 0.1.2.9000

## Major changes

* added inst/CITATION file    
* aligned ordering and names of arguments for both `boot_pvl` and `scan_pvl`    
* added examples for `boot_pvl`    
* added literature references to both `boot_pvl` and `scan_pvl`    


## Bug fixes

* corrected typo in vignette  



# qtl2pleio 0.1.2

## Major changes

* added tests & examples for `scan_pvl`  
* changed output of `scan_pvl` to a tibble  
* added `boot_n` function for use in bootstrap analyses  
* started using covariates in calls to `calc_covs`. Note that we still don't use genetic data when calling `calc_covs`.  
* deprecated `calc_loglik_bvlmm`



# qtl2pleio 0.1.1

## Major changes

* restructured `scan_pvl` to allow for more than two phenotypes. Now output is a dataframe.

# qtl2pleio 0.1.0

## Major changes

* Added a `NEWS.md` file to track changes to the package.
