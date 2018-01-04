#' Read boot output files - which contain lrt statistics - and assemble into a single table
#'
#' @param DIR directory path for directory that contains boot output files
#' @export

make_results_table <- function(DIR = "results-chtc/boot400-run1"){
  i <- 1
  out <- list()
  dir(DIR) -> fns
  for (fn in fns){
    read.table(file.path(DIR, fn)) -> foo
    dplyr::mutate(foo, filename = fn) -> foo2
    foo2 -> out[[i]]
    i <- i + 1
  }
  do.call("rbind", out) -> mytab
  return(mytab)
}
