# qtl2pleio 0.1.2.9001

## Major changes

* added a vignette for HTCondor & bootstrap analysis
* updated README.Rmd and README.md


# qtl2pleio 0.1.2.9000

## Major changes

* added inst/CITATION file    
* aligned ordering and names of arguments for both `boot_pvl` and `scan_pvl`    
* added examples for `boot_pvl`    
* added literature references to both `boot_pvl` and `scan_pvl`    




## Bug fixes

* corrected typo in vignette  



# qtl2pleio 0.1.2

* added tests & examples for `scan_pvl`  
* changed output of `scan_pvl` to a tibble  
* added `boot_n` function for use in bootstrap analyses  
* started using covariates in calls to `calc_covs`. Note that we still don't use genetic data when calling `calc_covs`.  
* deprecated `calc_loglik_bvlmm`



# qtl2pleio 0.1.1

* restructured `scan_pvl` to allow for more than two phenotypes. Now output is a dataframe.

# qtl2pleio 0.1.0

* Added a `NEWS.md` file to track changes to the package.
