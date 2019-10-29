## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.5.3
* win-builder (devel and release)

## R CMD check results  

0 errors | 0 warnings | 2 notes

* This is a new release.

* The package qtl2 is available through a non-CRAN repository.

## Replies to previous CRAN concerns

* I shortened the title in the file DESCRIPTION

* I stopped exporting some functions that are not needed directly by the user.


* This is a resubmission of an initially rejected package. I mistakenly included code that downloaded files to a new directory. I've since deleted those pieces of code.

* I omitted the abbreviation QTL in the DESCRIPTION file.


