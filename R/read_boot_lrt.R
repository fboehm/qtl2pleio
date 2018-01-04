#' Read LRT statistics for all of the scans (ie, 400 typically) that correspond to a single simulated trait file
#'
#' @param job_num trait file id number. An integer in 0:399 typically
#' @param total_number_of_scans number of scans ie number of 2d scans - for each trait file
#' @param nlrt_per_file positive integer specifying number of lrt statistics per boot output file
#' @param boot_run_num run number for the boot files. Used to create the directory path and file names
#' @export
read_boot_lrt <- function(job_num = 0, total_number_of_scans = 400,
                          nlrt_per_file = 10,
                          boot_run_num = 1
){
  stopifnot(job_num < total_number_of_scans)
  DIR <- paste0("~/Box Sync/attie/research-notebook/results-chtc/boot400-run", boot_run_num)
  nfiles <- total_number_of_scans / nlrt_per_file
  out <- list()
  job_nums <- job_num + (0:(nfiles - 1)) * total_number_of_scans
  fns <- paste0(job_nums, "_sim1_sim2-run", boot_run_num, ".txt")

  for (i in 1:nfiles){
    fn <- file.path(DIR, fns[i])
    read.table(fn) -> foo
    dplyr::mutate(foo, filename = fn) -> out[[i]]
  }
  do.call("rbind", out) -> mytab
  return(mytab)
}
