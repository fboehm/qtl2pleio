#' Assemble tibble from matrix of log-likelihood values
#'
#' @family profile log-likelihood tibble functions
#' @param loglik_mat a square matrix of log-likelihood values with rownames and colnames matching
#' @return tibble with 3 columns: marker1, marker2 and log-likelihood
#' @examples
#' llmat <- matrix(nrow = 3, ncol = 3, data = rgamma(9, 5))
#' rownames(llmat) <- paste0('m', 1:3)
#' colnames(llmat) <- paste0('m', 1:3)
#' transform_loglik_mat(llmat)
#' @export
transform_loglik_mat <- function(loglik_mat) {
    marker1 <- rep(rownames(loglik_mat), times = ncol(loglik_mat))
    marker2 <- rep(colnames(loglik_mat), each = nrow(loglik_mat))
    ll <- as.vector(loglik_mat)
    dat <- tibble::tibble(marker1, marker2, ll)
    return(dat)
}

#' Add physical map contents to tibble
#'
#' @family profile log-likelihood tibble functions
#' @param tib a tibble with 3 columns: marker1 name, marker2 name and log-likelihood values
#' @param pmap a physical map for a single chromosome
#' @return a tibble with 5 columns: marker1, marker2, log-likelihood, marker1_position, marker2_position
#' @examples
#' pm <- 1:3
#' names(pm) <- paste0('m', 1:3)
#' expand.grid(paste0('m', 1:3), paste0('m', 1:3)) -> foo
#' names(foo) <- c('marker1', 'marker2')
#' foo$ll <- rgamma(9, 5)
#' tibble::as_tibble(foo) -> tib
#' add_pmap(tib, pm)
#' @export
#' @importFrom rlang .data

add_pmap <- function(tib, pmap) {
    pmap_tib <- tibble::tibble(marker = names(pmap), marker_position = pmap)
    # now join
    tib2 <- dplyr::left_join(tib, pmap_tib, by = c(marker1 = "marker"))
    # rename(marker1_position = marker_position
    tib3 <- dplyr::rename(tib2, marker1_position = .data$marker_position)
    tib4 <- dplyr::left_join(tib3, pmap_tib, by = c(marker2 = "marker"))
    tib5 <- dplyr::rename(tib4, marker2_position = .data$marker_position)
    return(tib5)
}

#' Assemble a profile ll tibble
#'
#' @family profile log-likelihood tibble functions
#' @param tib a tibble, derived from loglik_mat and with pmap info
#' @param trace character vector to identify which profile log likelihood to calculate. Takes values 'profile1' or 'profile2'
#' @export
#' @return a tibble with 3 columns: marker position, profile log-likelihood, and trace (profile1 or profile2)
#'
#' @importFrom rlang .data

assemble_profile_tib <- function(tib, trace = "profile1") {
    if (trace == "profile1") 
        tib2 <- dplyr::group_by(tib, .data$marker1)
    if (trace == "profile2") 
        tib2 <- dplyr::group_by(tib, .data$marker2)
    tib3 <- dplyr::filter(tib2, .data$ll == max(.data$ll))
    tib4 <- dplyr::ungroup(tib3)
    if (trace == "profile1") {
        tib4$trace <- "profile1"
        tib4$marker_position <- tib4$marker1_position
    }
    if (trace == "profile2") {
        tib4$trace <- "profile2"
        tib4$marker_position <- tib4$marker2_position
    }
    tib5 <- dplyr::select(tib4, .data$marker_position, .data$ll, .data$trace)
    return(tib5)
}

#' Tidy the data frame outputted by scan_pvl for further analysis & plotting
#'
#' @family profile log-likelihood tibble functions
#' @param mytib outputted dataframe from scan_pvl
#' @param pmap physical map (in Mb) for exactly one chromosome, pmap$`5`, for example
#' @export
#' @importFrom rlang .data

tidy_scan_pvl <- function(mytib, pmap) {
    mytib <- tibble::as_tibble(mytib)
    mytib <- dplyr::rename(mytib, marker1 = .data$Var1, marker2 = .data$Var2, ll = .data$loglik)
    dat <- add_pmap(mytib, pmap)
    pl <- dplyr::filter(dat, .data$marker1 == .data$marker2)
    # pleio_ll <- dplyr::rename(pleio_ll, marker_position = .data$marker1_position)
    pl2 <- dplyr::rename(pl, marker_position = .data$marker1_position)
    # pleio_ll <- dplyr::mutate(pleio_ll, trace = 'pleio')
    pl2$trace <- "pleio"
    pleio_ll <- dplyr::select(pl2, .data$marker_position, .data$ll, .data$trace)
    # assemble pro1_ll
    pro1 <- assemble_profile_tib(dat, "profile1")
    pro2 <- assemble_profile_tib(dat, "profile2")
    # bind 3 tibbles
    foo <- dplyr::bind_rows(pleio_ll, pro1, pro2)
    foo$lod <- (foo$ll - max(pleio_ll$ll))/log(10)  # convert from base e to base 10
    dat <- dplyr::select(foo, .data$marker_position, .data$lod, .data$trace)
    return(dat)
}

#' Add intercepts to tidied loglikelihood tibble
#'
#' @family profile log-likelihood tibble functions
#' @param dat a tibble that results from tidy_scan_pvl acting on a log likelihood matrix
#' @param intercepts_univariate a vector of length 2 that contains the x coordinate values for the ordered univariate peaks
#' @export
#' @importFrom rlang .data

add_intercepts <- function(dat, intercepts_univariate) {
    dat$intercept_uni <- NA
    dat$intercept_uni[dat$trace == "profile1"] <- intercepts_univariate[1]
    dat$intercept_uni[dat$trace == "profile2"] <- intercepts_univariate[2]
    # determine pleiotropy peak location
    pleio_lod <- dplyr::filter(dat, .data$trace == "pleio")
    pleio_peak <- pleio_lod$marker_position[which.max(pleio_lod$lod)]
    dat$intercept_pleio <- NA
    dat$intercept_pleio[dat$trace == "pleio"] <- pleio_peak
    # determine bivariate peak locations
    dat$intercept_biv <- NA
    profile1_lod <- dplyr::filter(dat, .data$trace == "profile1")
    profile1_peak <- profile1_lod$marker_position[which.max(profile1_lod$lod)]
    profile2_lod <- dplyr::filter(dat, .data$trace == "profile2")
    profile2_peak <- profile2_lod$marker_position[which.max(profile2_lod$lod)]
    dat$intercept_biv[dat$trace == "profile1"] <- profile1_peak
    dat$intercept_biv[dat$trace == "profile2"] <- profile2_peak
    return(dat)
}

