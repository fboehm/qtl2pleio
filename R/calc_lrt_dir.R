#' Calculate likelihood ratio tests for all files in a single directory. Each file contains output of a two-dimensional scan.
#'
#' @param directory path to directory that contains the collection of files
#' @export

calc_lrt_dir <- function(directory){
  dir(directory) -> fns
  lrt <- numeric()
  i <- 1
  for (fn in fns){
    read.table(file.path(directory, fn)) -> scan_out
    calc_lrt(as.matrix(scan_out)) -> lrt[i]
    i <- i + 1
  }
  return(lrt)
}
