#' Plot tidied results of a pvl scan
#'
#' @param dat a tibble that is the tidied output of a pvl scan
#' @param phenames vector of names of the two ordered phenotypes
#' @param palette a vector containing strings for colors
#' @export
#' @importFrom rlang .data

plot_pvl <- function(dat, phenames, size = 3, shape = 17, palette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")){
  # find the

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
    ggplot2::geom_point(ggplot2::aes(x = dat$intercept,
                            y = c(rep(min(dat$lod), nrow(dat)))),
                        shape = shape,
                        size = size
                        ) +
    ggplot2::geom_point(ggplot2::aes(x = ))
}
