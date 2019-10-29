#' Plot tidied results of a pvl scan
#'
#' @param dat a profile lod tibble
#' @param units a character vector of length one to indicate units for physical or genetic map
#' @param palette a character vector of length 3 containing strings for colors
#' @param linetype a character vector of length 3 specifying the linetype values for the 3 traces
#' @export
#' @return a ggplot object with profile LODs

plot_pvl <- function(dat, units = "Mb",
                     palette = c("#999999", "#E69F00", "#56B4E9"),
                     linetype = c("solid", "longdash", "dotted"))
{
  dat %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(y = lod, x = marker_position,
                           color = trace, linetype = trace)) +
    ggplot2::scale_colour_manual(values = c("pleiotropy" = palette[1],
                                            "profile1" = palette[2],
                                            "profile2" = palette[3]),
                                 labels = unique(dat$Trace)) +
    ggplot2::scale_linetype_manual(values = c("pleiotropy" = linetype[1],
                                              "profile1" = linetype[2],
                                              "profile2" = linetype[3]),
                                   labels = unique(dat$Trace)) +
    ggplot2::labs(x = paste0("Chromosome position (", units, ")"), y = "LOD")
}

