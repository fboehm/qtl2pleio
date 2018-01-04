#' Create a pleiotropy peaks file with the marker indices - that designate the pleiotropy peak position - for each of a collection of trait files
#'
#' @param run_num run number, a positive integer
#' @param map a map for the corresponding chromosome
#' @param DIR a directory path for finding the 2d scan results
#' @export

write_pleio_peaks_file <- function(run_num, map = pm,
                                   DIR = paste0("results-chtc/pvl400-run", run_num)
){
  i <- 1
  peak_positions <- numeric()
  fns <- dir(DIR)
  for (fn in fns){
    as.matrix(read.table(file.path(DIR, fn))) -> out
    which.max(diag(out)) -> peak_positions[i]
    i <- i + 1
  }
  names(peak_positions) <- rownames(out)[peak_positions]
  rsid <- names(peak_positions)
  # load pmap, which is needed to find index of pleio peak.
  pleio_indices <- match(rsid, names(map))
  foo <- stringr::str_split(fns, pattern = "_")
  bar <- sapply(FUN = function(x)x[2], X = foo)
  foo2 <- stringr::str_split(bar, pattern = ".txt")
  jobnum <- as.numeric(sapply(FUN = function(x)x[1], X = foo2))
  # Look at the Ysim trait files to get their names
  Yfns <- paste0("Ysim-run", run_num, "_", jobnum, ".txt")
  tibble::tibble(pleio_indices, rsid, fns, jobnum, Yfns) -> tib
  # save tib
  write_csv(x = tib, paste0("pleio-peak-indices-table-run", run_num, ".csv"))
}
