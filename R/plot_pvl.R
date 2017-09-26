#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @param x1 x coordinate where trait 1 has its univariate peak
#' @param x2 x coordinate where trait 2 has its univariate peak
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat){
    # plot
    ggplot2::ggplot(dat, ggplot2::aes(y = dat$lod, x = dat$marker_position, colour = dat$trace)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "LOD", x = "marker position") +
    ggplot2::geom_vline(xintercept = x1) +
    ggplot2::geom_vline(xintercept = x2)
}
