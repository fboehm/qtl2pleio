#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat){
    # organize pleiotropy lod scores
    pleio_ll <- dplyr::filter(dat, .data$marker1 == .data$marker2)
    #pleio_ll <- dplyr::rename(pleio_ll, marker_position = .data$marker1_position)
    pleio_ll$marker_position <- pleio_ll$marker1_position
    #pleio_ll <- dplyr::mutate(pleio_ll, trace = "pleio")
    pleio_ll$trace <- "pleio"
    pleio_ll <- dplyr::select(pleio_ll, .data$marker_position, .data$ll, .data$trace)
    # assemble pro1_ll
    pro1_ll <- dplyr::group_by(dat, .data$marker1_position)
    pro1_ll <- dplyr::select(pro1_ll, .data$marker1, .data$ll, .data$marker1_position)
    #pro1_ll <- dplyr::rename(pro1_ll, marker = .data$marker1, marker_position = .data$marker1_position)
    pro1_ll$marker <- pro1_ll$marker1
    pro1_ll$marker_position <- pro1_ll$marker1_position

    #pro1_ll <- dplyr::mutate(pro1_ll, trace = "profile1")
    pro1_ll$trace <- "profile1"
    pro1_ll <- dplyr::filter(pro1_ll, .data$ll == max(.data$ll))
    pro1_ll <- dplyr::select(pro1_ll, .data$marker_position, .data$ll, .data$trace)
    # assemble pro2_ll
    pro2_ll <- dplyr::group_by(dat, .data$marker2_position)
    pro2_ll <- dplyr::select(pro2_ll, .data$marker2, .data$ll, .data$marker2_position)
    #pro2_ll <- dplyr::rename(pro2_ll, marker = .data$marker2, marker_position = .data$marker2_position)
    pro2_ll$marker <- pro2_ll$marker2
    pro2_ll$marker_position <- pro2_ll$marker2_position
    #pro2_ll <- dplyr::mutate(pro2_ll, trace = "profile2")
    pro2_ll$trace <- "profile2"
    pro2_ll <- dplyr::filter(pro2_ll, .data$ll == max(.data$ll))
    pro2_ll <- dplyr::select(pro2_ll, .data$marker_position, .data$ll, .data$trace)
    # bind 3 tibbles
    dat <- dplyr::bind_rows(pleio_ll, pro1_ll, pro2_ll)
    # plot
    ggplot2::ggplot(dat, ggplot2::aes(y = .data$ll - max(pleio_ll$ll), x = .data$marker_position, colour = .data$trace)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Log likelihood difference")
}
