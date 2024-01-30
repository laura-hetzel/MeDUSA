# *** Quality Metrics -----------------------------------------------------
#' MZ-OBJ Quality Metrics
#'
#' Get Statistical Mz Data\cr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#'
#' Dependencies : dplyr
#' @return Returns a dataframe of samples by metrics
#' @export
mz_quality_metrics <- function(input_mz_obj, cores = 2){
  no_mz <- dplyr::select(input_mz_obj, -mz)
  mz_metrics <- data.frame( name = colnames(no_mz),
                            median_mz = numeric(length(no_mz)),
                            min_mz = numeric(length(no_mz)),
                            max_mz = numeric(length(no_mz)),
                            n_peaks = numeric(length(no_mz)),
                            peaks_1k = numeric(length(no_mz)),
                            peaks_10k = numeric(length(no_mz)),
                            peaks_100k = numeric(length(no_mz)))

  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    tmp <- pbapply::pbapply(no_mz, 2,cl=cl,function(pop_sum){
      test <- input_mz_obj$mz[pop_sum>0]
      mz_metrics[colnames(pop_sum),] <-  c(median_mz = median(test),
                                           min_mz = min(test),
                                           max_mz = max(test),
                                           n_peaks = length(test),
                                           peaks_1k = sum(pop_sum > 1000),
                                           peaks_10k = sum(pop_sum > 10000),
                                           peaks_100k = sum(pop_sum >100000))
    })
    mz_metrics[-1] <- t(tmp)
    return(data.frame(mz_metrics))
  },
  finally={
    local.kill_threads(cl)
  })
}


# *** Metrics Quality Plot All -----------------------------------------------------
#' MZ_METRICS plot all statistics
#'
#' @param mz_metrics \cr
#'   DataFrame : from "sumR::mz_quality_metrics"
#' @export
mzmetrics_quality_plot_all <- function(mz_metrics){
  for (f in colnames(mz_metrics[-1])){
    mzmetrics_quality_plot(mz_metrics, f)
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
#' Dependencies : ggplot2, ggpubr, dplyr, parallel
#' @export
mzmetrics_quality_plot <- function(mz_metrics, focus, title = F, plot_dim = c(8,8)){
  if(title == F ){
    title <- focus
  }
  
  ggpubr::ggbarplot(mz_metrics, x = "name", y = focus,
                    title = title,
                    fill  = focus,
                    color = focus,    
                    xlab = "Sample_Name",
                    ylab = title) + 
    ggpubr::rotate(ylim = c( min(mz_metrics[[focus]][is.finite(mz_metrics[[focus]])]), 
                             max(mz_metrics[[focus]][is.finite(mz_metrics[[focus]])]) )) +
    ggplot2::scale_x_discrete(
      label = sapply(strsplit(mz_metrics$name, '_'), function(x) paste(x[-2:-1], collapse = '.')), 
      guide = ggplot2::guide_axis( n.dodge=3 )) +
    ggplot2::theme(legend.position = "none", axis.text.y=ggplot2::element_text(size=ggplot2::rel(0.5)))
  
  local.save_plot(paste("QualityCheck", title, local.mz_polarity_guesser(mz_metrics),sep="-"), dim = plot_dim)
}

# *** Metrics Quality Metadata Check -----------------------------------------------------
#' Does the metadata seem ok
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param meta \cr
#'   DataFrame : Input MZ-Obj
#'
#' Dependencies : dplyr
#' @export
mz_quality_meta_check <- function(input_mz_obj, meta){
  meta <- local.meta_polarity_fixer(input_mz_obj, meta)
  .subset_check <- function( listA, listB, msg ){
    missing <- listA[!(listA %in% listB)]
    if( length(missing) > 0 ){
      error <- append(error,paste0(msg ,missing))
    }
    error
  }

  error <- list()

  #Many more are expected for good science, but these are the ones that absence will actively break
  meta_columns_expected <- c("measurement",
                             "sample_name",
                             "type",
                             "phenotype",
                             "filtered_out")

  error <- .subset_check(meta_columns_expected,colnames(meta),"MetaColumns missing: ")
  error <- .subset_check(meta$sample_name,colnames(input_mz_obj),"MetaSamples not in data: ")
  error <- .subset_check(colnames(dplyr::select(input_mz_obj,-mz)), meta$sample_name,"DataSamples not in meta: ")

  sample_diff <- nrow(meta) - ncol(dplyr::select(input_mz_obj, -mz))
  if ( sample_diff != 0){
    error <- append(paste0(sample_diff,": more meta samples than data ones."),error)
  }


  if (length(error) > 0){
    stop(paste0("ERROR: sumR::mz_quality_meta_check: ", error, "\n"))
  }
}

# *** Quality Magic -----------------------------------------------------
#' MZ-OBJ Quality Magic
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param meta \cr
#'   DataFrame : Input MZ-Obj
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#'
#' Dependencies : dplyr
#' @return Returns a dataframe of samples by metrics
#' @export
mz_quality_magic <- function(input_mz_obj, meta, cores = 2){
  mz_quality_meta_check(input_mz_obj, meta)
  qc <- mz_quality_metrics(input_mz_obj, cores)
  mzmetrics_quality_plot_all(qc)
  qc
}
