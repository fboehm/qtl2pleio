#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @export

plot_pvl <- function(dat){
    # organize pleiotropy lod scores
    pleio_ll <- dplyr::filter(dat, marker1 == marker2)
    pleio_ll <- dplyr::rename(pleio_ll, marker_position = marker1_position)
    pleio_ll <- dplyr::mutate(pleio_ll, trace = "pleio")
    pleio_ll <- dplyr::select(pleio_ll, marker_position, ll, trace)
    # assemble pro1_ll
    pro1_ll <- dplyr::group_by(dat, marker1_position)
    pro1_ll <- dplyr::select(pro1_ll, marker1, ll, marker1_position)
    pro1_ll <- dplyr::rename(pro1_ll, marker = marker1, marker_position = marker1_position)
    pro1_ll <- dplyr::mutate(pro1_ll, trace = "profile1")
    pro1_ll <- dplyr::filter(pro1_ll, ll == max(ll))
    pro1_ll <- dplyr::select(pro1_ll, marker_position, ll, trace)
    # assemble pro2_ll
    pro2_ll <- dplyr::group_by(dat, marker2_position)
    pro2_ll <- dplyr::select(marker2, ll, marker2_position)
    pro2_ll <- dplyr::rename(marker = marker2, marker_position = marker2_position)
    pro2_ll <- dplyr::mutate(trace = "profile2")
    pro2_ll <- dplyr::filter(ll == max(ll))
    pro2_ll <- dplyr::select(marker_position, ll, trace)
    # bind 3 tibbles
    dat <- dplyr::bind_rows(pleio_ll, pro1_ll, pro2_ll)
    # plot
    ggplot2::ggplot(dat, aes(y = ll - max(pleio_ll$ll), x = marker_position, colour = trace)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Log likelihood difference")
}
