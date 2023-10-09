# *** Quality Metrics -----------------------------------------------------
#' MZ-OBJ Quality Metrics
#'
#' Get Statistical Mz Data\cr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#'
#' Dependencies : dplyr, parallel
#' @return Returns a dataframe of samples by metrics
#' @export
mz_quality_metrics <- function(input_mz_obj){
  no_mz <- dplyr::select(input_mz_obj, -mz)
  mz_metrics <- data.frame( name = colnames(no_mz),
                            median_mz = numeric(length(no_mz)),
                            min_mz = numeric(length(no_mz)),
                            max_mz = numeric(length(no_mz)),
                            n_peaks = numeric(length(no_mz)),
                            peaks_1k = numeric(length(no_mz)),
                            peaks_10k = numeric(length(no_mz)),
                            peaks_100k = numeric(length(no_mz)))

  cl <- local.export_thread_env(cores, deparse(sys.calls()[[sys.nframe()]]))
  tryCatch({
    tmp <- pbapply::pbapply(no_mz, 2,cl=cl,function(pop_sum){
      test <- input_mz_obj$mz[pop_sum>0]
      mz_metrics[colnames(pop_sum),] <-  c(median_mz = median(test),
                                           min_mz = min(test),
                                           maz_mz = max(test),
                                           n_peaks = length(test),
                                           peaks_1k = sum(pop_sum > 1000),
                                           peaks_10k = sum(pop_sum > 10000),
                                           peaks_100k = sum(pop_sum >100000))
  })
    return(mz_metrics[-1] <- t(tmp))
  },
  finally={
    if (cores > 1 || !is.null(cl)) {
      parallel::stopCluster(cl)
      showConnections()
    }
  })
}


# *** Metrics Quality Plot All -----------------------------------------------------
#' MZ_METRICS plot all statistics
#'
#' @param mz_metrics \cr
#'   DataFrame : from "sumR::mz_quality_metrics"
#' @export
mz_metrics_quality_plot_all <- function(mz_metrics){
  for (f in colnames(mz_metrics[-1]){
    mz_metrics_quality_plot(mz_metrics, f)
    local.save_plot(paste("Quality Check ",f,local.mz_polarity_guesser(input_mz_obj),sep="-"))
  }
}

# *** Metrics Quality Plot -----------------------------------------------------
#' MZ_METRICS plot a single statistic
#'
#' @param mz_metrics \cr
#'   DataFrame : from "sumR::mz_quality_metrics"
#' @param focus \cr
#'   String : which metric to plot (i.e. median_mz)
#' @param title \cr
#'   String : title of plot
#'
#' Dependencies : ggplot2, dplyr, parallel
#' @export
mz_metrics_quality_plot <- function(mz_metrics, focus, title = focus){
  ggplot() +
    geom_line(data = mz_metrics,
              aes_string(x = "name" , y = focus, group = "1")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(title)
}
