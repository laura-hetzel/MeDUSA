# *** Quality Metrics -----------------------------------------------------
#' Run basic statictics over an mz_obj
#'
#' @description
#' The methods of single cell mass spectrometry have multiple challenges
#' associated with them, some of which may lead to a poor measurement that
#' does not accurately reflect the biological status of the sample. It is
#' recommended that these samples are removed from the data set. Before removal
#' from the data set, the entire data set can undergo a quality check. The user is
#' expected to be able to identify outliers from the data provided with the
#' mz_quality_metrics function. It is recommended that the user check the raw
#' files of identified outliers and confirm a likely poor measurement before
#' removing the file from the data set.
#'
#' MZ-OBJ Quality Metrics
#'
#' Get Statistical Mz Data\cr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param cores
#'   Integer: Can I has multithreading? (Requires parallel)
#'
#' Dependencies : dplyr
#' @examples
#'   mz_quality_metrics(input_mz_obj) : Perform quality metrics
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
#'   DataFrame : from "MeDUSA::mz_quality_metrics"
#' @examples
#'   mzmetrics_quality_plot_all(mz_metrics) : Plot all metrics
#' @return plots of all quality metrics
#' @export
mzmetrics_quality_plot_all <- function(mz_metrics){
  for (f in colnames(mz_metrics[-1])){
    mzmetrics_quality_plot(mz_metrics, f)
  }
}

# *** Metrics Quality Plot -----------------------------------------------------
#' MZ_METRICS plot a single statistic
#'
#' @description
#' Plot a single quality statistic;
#' @param mz_metrics \cr
#'   DataFrame : from "MeDUSA::mz_quality_metrics"
#' @param focus \cr
#'   String : which metric to plot:
#'     ( median_mz min_mz max_mz n_peaks peaks_1k peaks_10k peaks_100k )
#' @param title \cr
#'   String : title of plot; defaults to set to "focus"
#' @param plot_dim \cr
#'   c(Int,Int) : Dimensions of plot
#' Dependencies : ggplot2, ggpubr, dplyr, parallel
#' @examples
#'   mzmetrics_quality_plot(mz_metrics, 'peaks_1k') : Plot all metrics
#'   mzmetrics_quality_plot(mz_metrics, 'n_peaks', title = "Custom_npeaks", plot_dim = c(10,10)) 
#'     : Plot n_peaks with a custom title and dimensions
#' @return plots of a provided quality metrics
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
#' @description
#' The metadata file is used to filter the mz_obj to easily idetify samples by phenotype,
#' type (i.e. blank/cell/media), fltered_out (blacklisting), or anything else that the user
#' chooses to differentiate samples by.
#' Metadata requirements are:
#'   Includes at least : data.frame(sample_name("character),
#'                                  type("character"),
#'                                  phenotype("character"),
#'                                  filtered_out("logical"))
#'   Ignoring polarity in mz_obj. All mz_obj columns must be found in metadata$sample_name and vice-versa.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param meta \cr
#'   DataFrame : metadata object
#'   Dependencies : dplyr
#' @examples
#'   mz_quality_meta_check(input_mz_obj, meta) : run meta checks
#' @return
#'   null (only errors if metadata mismatches sample data)
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

  #Many more can be expected for "good science", but these are the ones that absence will cause issue
  meta_columns_expected <- c( "sample_name",
                              "type",
                              "phenotype",
                              "filtered_out" )


  error <- .subset_check(meta_columns_expected,colnames(meta),"MetaColumns missing: ")
  error <- .subset_check(meta$sample_name,colnames(input_mz_obj),"MetaSamples not in data: ")
  error <- .subset_check(colnames(dplyr::select(input_mz_obj,-mz)), meta$sample_name,"DataSamples not in meta: ")
  if ( class(meta$filtered_out) != 'logical' ){
    error <- append("metadata$filtered_out column must be a logical (boolean)")
  }

  sample_diff <- nrow(meta) - ncol(dplyr::select(input_mz_obj, -mz))
  if ( sample_diff != 0){
    error <- append(paste0(sample_diff,": more meta samples than data ones."),error)
  }

  if (length(error) > 0){
    stop(paste0("ERROR: MeDUSA::mz_quality_meta_check: ", error, "\n"))
  }
}

# *** Quality Magic -----------------------------------------------------
#' MZ-OBJ Quality Magic
#' Run all the quality things. Quality metrics, meta checks & generate plots
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param meta \cr
#'   DataFrame : metadata object
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#'
#' Dependencies : dplyr
#' @examples
#'   mz_quality_magic(input_mz_obj, meta) : run all the things
#' @return Returns a dataframe of samples by metrics
#' @export
mz_quality_magic <- function(input_mz_obj, meta, cores = 2){
  mz_quality_meta_check(input_mz_obj, meta)
  qc <- mz_quality_metrics(input_mz_obj, cores)
  mzmetrics_quality_plot_all(qc)
  qc
}
