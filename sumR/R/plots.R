#' Volcano plot
#'
#' @param data 
#' @param xvalues 
#' @param yvalues 
#' @param title 
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
volcanoPlot <- function(data, xvalues, yvalues, title){
  ggplot(data) +
    geom_point(aes(x=xvalues, y=-log10(yvalues), colour=significant)) + ## color by significant of fdr <0.1
    ggtitle(title) +
    xlab("log2 fold change") +
    ylab("-log10 nominal p.value") +
    #scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "right",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    scale_color_aaas() +
    theme_pubr() +
    labs_pubr()
}