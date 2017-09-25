#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat){
    # plot
    ggplot2::ggplot(dat, ggplot2::aes(y = dat$lod, x = dat$marker_position, colour = dat$trace)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "LOD", x = "marker position")
}
