#' Tidy the matrix of log likelihood values for further analysis & plotting
#'
#' @param loglik_mat a (square) matrix of log likelihood values
#' @param pmap physical map (in Mb) for exactly one chromosome, pmap$`5`, for example
#' @export
#' @importFrom rlang .data

tidy_scan_pvl <- function(loglik_mat, pmap){
  # Assumes that we have rownames and columns names assigned to loglik_mat
  marker1 <- rep(rownames(loglik_mat), times = ncol(loglik_mat))
  marker2 <- rep(colnames(loglik_mat), each = nrow(loglik_mat))
  ll <- as.vector(loglik_mat)
  dat <- tibble::tibble(marker1, marker2, ll)
  pmap_tib <- tibble::tibble(marker = names(pmap), marker_position = pmap)
  # now join
  foo <- dplyr::left_join(dat, pmap_tib, by = c("marker1" = "marker"))
  #  rename(marker1_position = marker_position
  foo$marker1_position <- foo$marker_position
  foo <- foo[, -4]
  foo <- dplyr::left_join(foo, pmap_tib, by = c("marker2" = "marker"))
  #  rename(marker2_position = marker_position)
  foo$marker2_position <- foo$marker_position
  dat <- foo[, -5]
  pleio_ll <- dplyr::filter(dat, .data$marker1 == .data$marker2)
  #pleio_ll <- dplyr::rename(pleio_ll, marker_position = .data$marker1_position)
  pleio_ll$marker_position <- pleio_ll$marker1_position
  pleio_ll <- pleio_ll[ , -c(4, 5)]
  #pleio_ll <- dplyr::mutate(pleio_ll, trace = "pleio")
  pleio_ll$trace <- "pleio"
  pleio_ll <- dplyr::select(pleio_ll, .data$marker_position, .data$ll, .data$trace)
  # assemble pro1_ll
  pro1_ll <- dplyr::group_by(dat, .data$marker1_position)
  pro1_ll <- dplyr::filter(pro1_ll, .data$ll == max(.data$ll))
  pro1_ll <- dplyr::ungroup(pro1_ll)
  pro1_ll$trace <- "profile1"
  pro1_ll$marker_position <- pro1_ll$marker1_position
  pro1_ll <- dplyr::select(pro1_ll, .data$marker_position, .data$ll, .data$trace)
  # assemble pro2_ll
  pro2_ll <- dplyr::group_by(dat, .data$marker2_position)
  pro2_ll <- dplyr::filter(pro2_ll, .data$ll == max(.data$ll))
  pro2_ll <- dplyr::ungroup(pro2_ll)
  pro2_ll$trace <- "profile2"
  pro2_ll$marker_position <- pro2_ll$marker2_position
  pro2_ll <- dplyr::select(pro2_ll, .data$marker_position, .data$ll, .data$trace)
  # bind 3 tibbles
  foo <- dplyr::bind_rows(pleio_ll, pro1_ll, pro2_ll)
  foo$lod <- (foo$ll - max(pleio_ll$ll)) / log(10) # convert from base e to base 10
  dat <- dplyr::select(foo, .data$marker_position, .data$lod, .data$trace)
  # define pleiotropy intercept value
  #pleio_int <- pleio_ll$marker_position[which.max(pleio_ll$ll)]
  # define intercept column
  #dat$intercept <- NA
  #dat$intercept[dat$trace == "profile1"] <- intercepts[1]
  #dat$intercept[dat$trace == "profile2"] <- intercepts[2]
  #dat$intercept[dat$trace == "pleio"] <- pleio_int
  return(dat)
}

#' Add intercepts to tidied loglikelihood tibble
#'
#' @param dat a tibble that results from tidy_scan_pvl acting on a log likelihood matrix
#' @param intercepts_univariate a vector of length 2 that contains the x coordinate values for the ordered univariate peaks
#' @export
#' @importFrom rlang .data

add_intercepts <- function(dat, intercepts_univariate){
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

