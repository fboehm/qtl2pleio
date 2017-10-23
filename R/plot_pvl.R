#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @param phenames vector of names of the two ordered phenotypes
#' @param size number indicating the plot size for plot characters
#' @param shape_uni number indicating the pch for univariate peaks
#' @param shape_pleio number indicating the pch for pleiotropy peak
#' @param shape_biv number indicating the pch for bivariate peaks
#' @param palette a vector containing strings for colors
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat, phenames, size = 3, shape_uni = 17, shape_pleio = 23, shape_biv = 18, palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")){
  # plot
    ggplot2::ggplot(dat, ggplot2::aes(y = dat$lod, x = dat$marker_position, colour = dat$trace)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_manual(values = palette[1:3],
                                 breaks = c("pleio", "profile1", "profile2"),
                                 labels = c("Pleiotropy", phenames[1], phenames[2])
                                 ) +
    ggplot2::labs(y = "LOD", x = "Marker position") +
    #ggplot2::geom_vline(ggplot2::aes(xintercept = dat$intercept, colour = dat$trace)) +
    ggplot2::guides(colour =
               ggplot2::guide_legend(
                 title = "Trace",
                 title.theme = ggplot2::element_text(
                   size = 15,
                   face = "italic",
                   #colour = "red",
                   angle = 0
                 )
               )
    ) +
    ggplot2::geom_point(ggplot2::aes(x = dat$intercept_uni,
                            y = c(rep(min(dat$lod), nrow(dat)))),
                        shape = shape_uni,
                        size = size,
                        na.rm = TRUE
                        ) +
    ggplot2::geom_point(ggplot2::aes(x = dat$intercept_pleio,
                                     y = c(rep(min(dat$lod), nrow(dat)))),
                        shape = shape_pleio,
                        size = size,
                        na.rm = TRUE
                        ) +
    ggplot2::geom_point(ggplot2::aes(x = dat$intercept_biv,
                                     y = c(rep(min(dat$lod), nrow(dat)))),
                        shape = shape_biv,
                        size = size,
                        na.rm = TRUE
    )
}
