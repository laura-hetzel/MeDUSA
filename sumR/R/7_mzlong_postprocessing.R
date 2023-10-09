# *** Log Transform -----------------------------------------------------
#' MZLONG-OBJ Log Transform
#'
#' Convert intensities of MZLong to log2
#'  For normal "mz_obj" use log2(mz_obj)
#'  - Requires: ggplot2
#'
#' @param input_mzlong \cr
#'   DataFrame : Input MZLONG-Obj
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#' @returns mzLong_obj
#'
#' @export
mzlong_pp_log_transform <- function(input_mzlong, plot = TRUE){
  input_mzlong$log <- log2(input_mzlong$intensity)
  if(plot){
    ggplot(input_mzlong, aes(x = sample, y = log)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste("Normalized, Filtered,",
        local.mz_polarity_guesser(input_mzlong)
        , sep=""))
    local.save_plot(paste("LogTransform",local.mz_polarity_guesser(input_mzlong),sep="-"))
  }
  input_mzlong
}
