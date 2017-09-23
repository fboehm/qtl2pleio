#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat){
    # plot
    ggplot2::ggplot(dat, ggplot2::aes(y = .data$ll - max(pleio_ll$ll), x = .data$marker_position, colour = .data$trace)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Log likelihood difference")
}
