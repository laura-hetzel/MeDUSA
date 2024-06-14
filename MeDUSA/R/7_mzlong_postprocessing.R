# *** Log Transform -----------------------------------------------------
#' MZLONG-OBJ Log Transform
#'
#' The mzlong_post_log function log2 transforms the intensities of the data set,
#' after the data set has been restructured via the pivot_longer function. Log
#' transformation of the data is highly recommended so that the observed m/z
#' relationships are proportional and not additive, making the statistical
#' analysis and interpretation more biologically relevant. Note that this
#' function works with the pivot_longer object, not the standard mz_object.
#'
#'  - Requires: ggplot2
#'
#' @param input_mzlong \cr
#'   DataFrame : Input MZLONG-Obj
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#' @returns mzLong_obj
#'
#' @returns MzLong_obj
#'
#' @export
mzlong_post_log <- function(input_mzlong, plot = TRUE){
  input_mzlong$log <- log2(input_mzlong$intensity)
  if(plot){
    tmp <- select(input_mzlong, -intensity)
    tmp$intensity <- tmp$log
    mzlongpost.plot(tmp, local.mz_polarity_guesser(input_mzlong), "PivotLongLog")
  }
  input_mzlong
}

mzlongpost.plot <- function(input_mzlong, polarity, plot_title){
  ggpubr::ggscatter(
    input_mzlong,
    x      = "sample_name", y    = "intensity",
    xlab   = "Sample_Name", ylab = "Intensity",
    shape  = 20,            size = 1,
    color  = "intensity",
    title  = paste(plot_title, polarity, sep=": "),
    rotate = T ) +
    ggplot2::scale_x_discrete(
      label = sapply(strsplit(input_mzlong$sample_name, '_'), function(x) paste(x[-2:-1], collapse = '.')),
      guide = ggplot2::guide_axis( n.dodge=3)) +
    ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=ggplot2::rel(0.5)))

  local.save_plot(paste(plot_title ,polarity, sep="-"))

  ggpubr::ggboxplot(
    input_mzlong,
    x             = "sample_name", y            = "intensity",
    xlab          = "Sample_Name", ylab         = "Intensity",
    outlier.shape = 3,             outlier.size = 0.1,
    title         = paste(plot_title, "box:", polarity, sep=" "),
    legend        = "none",
    rotate        = T
    ) +
    ggplot2::scale_x_discrete(
      label = sapply(strsplit(input_mzlong$sample_name, '_'), function(x) paste(x[-2:-1], collapse = '.')),
      guide = ggplot2::guide_axis( n.dodge=3)) +
    ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=ggplot2::rel(0.5)))

  local.save_plot(paste(plot_title, "box" ,polarity, sep="-"))
}
