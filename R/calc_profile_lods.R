
#' Add physical map contents to tibble
#'
#' @family profile lod tibble functions
#' @param tib a tibble with 3 columns: marker, trace, and profile lod values, typically outputted by calc_profile_lods()
#' @param pmap a physical map for a single chromosome
#' @return a tibble with 4 columns: marker, trait, profile_lod, marker_position
#' @examples
#' pm <- 1:3
#' names(pm) <- as.character(paste0('m', 1:3))
#' expand.grid(paste0('m', 1:3), paste0('m', 1:3)) %>%
#'     tibble::as_tibble() %>%
#'     dplyr::mutate(log10lik = rgamma(9, 5)) %>%
#'     calc_profile_lods() %>%
#'     add_pmap(pm)
#' @export
#' @importFrom rlang .data

add_pmap <- function(tib, pmap) {
    tibble::tibble(marker = as.character(names(pmap)), marker_position = pmap) %>%
        dplyr::inner_join(tib)
}





#' Calculate profile lods for all traits
#'
#' @param scan_pvl_out tibble outputted from scan_pvl
#' @export
#' @return a tibble with 3 columns, indicating 'marker identity, trace (pleiotropy or profile1, profile2, etc.), and value of the profile lod (base 10) for that trace at that marker.

calc_profile_lods <- function(scan_pvl_out){
    nc <- ncol(scan_pvl_out)
    out <- list()
    for (i in 1:(nc - 1)) {
        out[[i]] <- scan_pvl_out %>%
            dplyr::group_by_at(i) %>%
            dplyr::summarise(profile = max(log10lik)) %>%
            dplyr::mutate(trait = paste0("tr", i)) %>%
            dplyr::rename(marker = paste0("Var", i))
    }
    # Calculate the "pleiotropy" trace
    pleio_trace <- scan_pvl_out %>%
        dplyr::filter_all(dplyr::all_vars(. == Var1 | is.numeric(.))) %>%
        dplyr::rename(profile = log10lik) %>%
        dplyr::select(Var1, profile) %>%
        dplyr::rename(marker = Var1) %>%
        dplyr::mutate(trait = "pleiotropy")
    #
    out %>%
        dplyr::bind_rows() %>%
        dplyr::bind_rows(pleio_trace) %>%
        dplyr::mutate(pleio_max = max(pleio_trace$profile)) %>%
        dplyr::mutate(profile_lod = profile - pleio_max) %>%
        dplyr::select(marker, trait, profile_lod)
}

