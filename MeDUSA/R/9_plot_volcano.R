# *** Volcano Plot -----------------------------------------------------
#' welch: create volcano plot
#'
#' @description
#' A volcano plot plots the fold change on the x-axis and the p-value of the
#' t-test on the y-axis. The function plot_volcano assumes the user has already
#' performed fold change and t-test, so these values must be passed into the
#' function.
#'
#' Compare two phenotypes
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param phenotype_a \cr
#'   DataFrame: Filtered metadata for phenotypeA
#' @param phenotype_b \cr
#'   DataFrame: Filtered metadata for phenotypeB
#'
#' @export
plot_volcano <- function(welch, fold_change, title = "Volcano_Plot"){
  if ( sum(round(welch$mz,6) != round(fold_change$mz,6)) > 0){
    stop("ERROR: plot_volcano : Welch mz does not match FoldChange mz")
  }
  df <- data.frame("mz"   =  welch$mz,
                   "p"    = -log10(welch$p),
                   "fold" =  fold_change$fold,
                   "diff" =  rep("NONE",nrow(welch)))

  df$diff[df$fold > 0.6 & df$p > -log10(0.05)] <- "UP"
  df$diff[df$fold < -0.6 & df$p > -log10(0.05)] <- "DOWN"

ggpubr::ggscatter( data = df,
                    x = "fold",y="p",
                    title = title,
                    shape = "diff",
                    color ="diff", palette = c("blue", "black", "red")) +
                    geom_vline(xintercept = c(-0.6, 0.6), col = "green") +
                    geom_hline(yintercept = -log10(0.05), col = "green")

  local.save_plot(title)
}
