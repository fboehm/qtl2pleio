#' Find jobs ids of missing files from CHTC
#'
#' @param DIR directory for results files
#' @param njobs number of jobs expected - typically 400 or 16000
#' @param outfile file path for writing missing job ids
#' @export
#'
find_missing_jobs <- function(DIR = "results-chtc/boot400-run1", njobs = 16000, outfile = "bad-jobs-boot400-run1.txt"){
  dir(DIR) -> fns
  stringr::str_split(fns, pattern = "_") -> foo
  sapply(FUN = function(x)x[1], foo) -> foo2
  stringr::str_split(foo2, pattern = ".txt") -> foo3
  unlist(foo3) -> foo4
  as.numeric(foo4) -> bar
  table(c(0:(njobs - 1), bar)) -> mytab
  as.numeric(names(mytab[mytab == 1])) -> badjobs
  write.table(badjobs, file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
