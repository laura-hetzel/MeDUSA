#' @title (LEGACY, use read_mzml instead) Read MS data
#' @param path path to the file
#' @return Data.frame with ordered mz values
#' @export
#' @importFrom stringr str_remove
#' @importFrom dplyr bind_rows
read_msdata <- function(path = "data") {
  files <- list.files(path = "data", full.names = T)
  file_names <- str_remove(string = files, pattern = ".txt")
  file_list <- lapply(setNames(files, file_names), function(x) {
    read.table(x, col.names = c("mz", "intensity"))
  })
  nn <- sapply(file_list, nrow) # each sample's row number
  non_empty <- nn != 0L # if there's any empty sample; return logic vector
  ms_data <- bind_rows(file_list[non_empty], .id = "sample")

  return(ms_data[order(ms_data$mz), ])
}


#' @title Read DIMS data from an mzML file
#' @param mzml Path to a mzml file. If omitted,
#' the package example will be used
#' @importFrom MSnbase readMSData
#' @examples
#' read_mzml()
#' @export
#' @return XCMS object with spectra
read_mzml <- function(mzml = NULL){
  if (is.null(mzml)) {
    dir <- system.file("extdata", package = "sumR")
    mzml <- list.files(dir, pattern = ".*negative.*.mzML$", full.names = T)
  }
  tryCatch({
    pd <- data.frame(sampleNames = mzml, class = "sample")
    blanks <- grep("Blank", mzml)
    pd$class[blanks] <- "blank"
    data <- MSnbase::readMSData(mzml, mode = "onDisk",
                                pdata = new("NAnnotatedDataFrame", pd))
    return(data)
  }, error = function(x){
    stop("Cannot parse the given file, not recognized as a proper .mzML file")
  })
}

#' @title Combine spectra and centroid
#' @param data XCMS object
#' @importFrom tools file_path_sans_ext
#' @importFrom xcms smooth pickPeaks
#' @importFrom utils getFromNamespace
#' @import MSnbase
#' @export
combine_spectra_centroid <- function(data){
  combineSpectra <- utils::getFromNamespace("combineSpectra", "MSnbase")
  list_data <- split(data, f = data@featureData@data$polarity)
  names(list_data) <- c("negative", "positive")
  sapply(names(list_data), function(polarity){
    mzml <- data@phenoData@data$sampleNames
    file <- sprintf("%s_centroided_%s.mzML",
                    basename(tools::file_path_sans_ext(mzml)),
                    polarity)
    if (file.exists(file)) {
      file.remove(file)
    }

    list_data[[polarity]] %>%
      combineSpectra(method = function(x){
        meanMzInts(x, intensityFun = combine_intensity)#base::max)
      }) %>%
      smooth() %>%
      pickPeaks() %>%
      writeMSData(file)
    file
  })
}


combine_intensity <- function(x){
  m <- mean(x[x > 0])
  ifelse(is.nan(m), 0, m)
}

#' @title DIMS pipeline
#' @importFrom xcms MSWParam findChromPeaks groupChromPeaks MzClustParam
#' fillChromPeaks FillChromPeaksParam chromPeaks
#' @importFrom pbapply pblapply
#' @importFrom logging loginfo
#' @examples
#' dims_pipeline(read_mzml())
#' @export
dims_pipeline <- function(data, plot_md = FALSE){
  suppressWarnings({
    loginfo("Running DIMS pipeline")
    groups <- seq_len(nrow(data@featureData))
    per_injection <- split(data, f = groups)
    loginfo(sprintf("Identifying peaks in %s spectra", length(groups)))

    params <- xcms::MSWParam(SNR.method = "data.mean", winSize.noise = 500)
    data <- do.call(c, pblapply(per_injection, cl = 1, function(inj){
      dat <-  suppressMessages(xcms::findChromPeaks(inj, param = params))
      # Filter using exclusing list
      df <- as.data.frame(xcms::chromPeaks(dat))
      df <- md_pipeline(df[,c(1,2,3)], plot_md)
      # filter the peaks using the md
      # xcms::filterChromPeaks(dat)
      dat
    }))


    clust_param <- xcms::MzClustParam(sampleGroups = groups)
    data <- xcms::groupChromPeaks(data, param = clust_param)
    data <- fillChromPeaks(data, param = FillChromPeaksParam())
    loginfo("DIMS pipeline complete")
  })
  data
}



#' @title Process features obtained from XCMS
#' @param data XCMS object with grouped peaks
#' @importFrom xcms featureValues featureDefinitions
#' @importFrom pmp filter_peaks_by_blank mv_imputation
#' @importFrom Biobase pData
#' @examples
#' data <- dims_pipeline(read_mzml())
#' feature_processing(data)
#' @export
feature_processing <- function(data){
  intensity_df <- featureValues(data)
  classes <- Biobase::pData(data)$class
  intensity_df <- filter_peaks_by_blank(intensity_df, fold_change = 5,
                        classes = classes,
                        blank_label = "blank",
                        remove_samples = TRUE)
  df_imp <- mv_imputation(df = intensity_df, method = "sv")

  mz <- xcms::featureDefinitions(data)[rownames(df_imp), ]$mzmed
  df <- cbind(mz = mz, df_imp)
  as.data.frame(df)
}


#' ppm calculation
#'
#' @description calculates the parts per million error between two different
#' masses
#' @export
ppm_calc <- function(mass1, mass2) {
  ((mass1 - mass2) / mass1) * 1e6
}

#' Align check
#'
#' @description align_check makes sure that all the m/z values are
#' aligned/binned correctly
#'
#' align_check takes 2 arguments (1 optional), a dataframe of peaks with mz
#' column, and an optinal argument for the coordinates of the plot to be zoomed
#' in on
#'
#' align_check outputs a list of three elements:
#' 1- Dataframe of 1 column containing the ppm error values
#' 2- boxplot of the ppm erro values with xcoords zoomed in to -20,0 (default)
#' 3- table of summary stats of ppm error values
#'
#' @examples
#' @export
#' @importFrom tibble as_tibble_col
#' @import ggplot2
align_check <- function(data_frame_fn, xcoords = c(-20, 0)) {
  odd_ind_fn <- seq(3, length(data_frame_fn$mz), 2)
  even_ind_fn <- seq(2, length(data_frame_fn$mz), 2)
  ppm_err_fn <- ppm_calc(data_frame_fn$mz[even_ind_fn],
                         data_frame_fn$mz[odd_ind_fn]) %>%
    as_tibble_col(column_name = "ppm_error")
  ppm_err_plot_fn <- ggplot(ppm_err_fn, aes(x = ppm_error)) +
    geom_boxplot() +
    coord_cartesian(xlim = xcoords) +
    ggtitle("ppm error boxplot") +
    theme_classic(base_size = 20)
  ppm_err_summary_fn <- summary(ppm_err_fn)
  ppm_err_list <- list(
    ppm_error_df = ppm_err_fn,
    boxplot = ppm_err_plot_fn,
    summary_stats = ppm_err_summary_fn
  )
  # paste(nrow(data_frame) - ) add a print
  return(ppm_err_list)
}

#' Title
#'
#' @param mass
#' @param samples
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
condition <- function(mass, samples, tolerance) {
  if (anyDuplicated(samples)) {
    return(NA)
  }
  mean_mass <- mean(mass)
  if (any(abs(mass - mean_mass) / mean_mass > tolerance)) {
    return(NA)
  }

  return(mean_mass)
}

#' Background removal
#'
#' @param dataframe
#' @param filter_type
#' @param blank_thresh
#' @param nsamples_thresh
#' @param blank_regx
#' @param filtered_df
#'
#' @return

#'
#' @examples
#' @importFrom dplyr select contains if_else filter
#' @export
blank_subtraction <- function(dataframe, filter_type = "median",
                              blank_thresh = 1, nsamples_thresh = 1,
                              blank_regx = "blank", filtered_df = FALSE) {
  ## Creates a dataframe with only blanks
  blanks <- select(dataframe, contains(blank_regx))
  ## Creates a dataframe without blanks
  samples <- select(dataframe, !contains(blank_regx))

  ## Adds column with median/max/etc. of each row times a specified
  ## threshold percentage
  blanks$threshold <- apply(blanks, 1, filter_type) * blank_thresh

  ## Calculates for each row how many data points don't exceed the initial
  ## threshold if that amount meets or exceeds a certain threshold percentage
  ## of the total amount of samples it will be marked as failed, if it remains
  ## below the threshold percentage it will be marked as passed. Finally, a
  ## column is added to the sample dataframe that shows which m/z values passed
  ## or failed
  for (i in seq_len(nrow(samples))) {
    below_thresh <- sum(samples[i, seq_len(ncol(samples))] <=
                          blanks$threshold[i])
    boundary <- below_thresh / ncol(samples) * 100 >= nsamples_thresh
    samples$index[i] <- if_else(boundary, "blank_fail", "blank_pass")
  }

  ## Dataframe that only displays the passed values, filtering out failed values
  filtered_samples <- filter(samples, index == "blank_pass")
  ## Prints the number of passed values
  print_pass <- nrow(filter(samples, index == "blank_pass"))
  ## Prints the number of failed values
  print_fail <- nrow(filter(samples, index == "blank_fail"))

  if (filtered_df == TRUE) {
    ## Makes the function return a list with a dataframe containing m/z
    ## rows that passed the threshold, the number of passed values and the
    ## number of failed values
    filtered_list <- list(
      "filtered_df" = filtered_samples,
      "n_passed" = print_pass,
      "n_failed" = print_fail
    )
    return(filtered_list)
  }

  ## Makes the function return a list with a dataframe containing m/z rows that
  ## both passed and failed to pass the threshold, the number of passed values
  ## and the number of failed values
  unfiltered_list <- list(
    "unfiltered_df" = samples,
    "n_passed" = print_pass,
    "n_failed" = print_fail
  )
  return(unfiltered_list)
}

#' @title Calculate nominal mass from a formula
#' @param formulas formulas in string format like 'c6h12o6'
#' @importFrom stringr str_extract_all
#' @importFrom PeriodicTable mass
calculate_nominal_mass <- function(formulas){
  sapply(formulas, function(formula){
    vec <- as.vector(str_extract_all(formula, "([A-Za-z][[:digit:]]+|[A-Z][a-z]{0,1})",
                                     simplify = T))
    sum(sapply(vec, function(comb){
      if (!grepl("[[:digit:]]", comb)) {
        atom <- comb
        mult <- 1
      } else {
        atom <- str_extract(comb, "[[:alpha:]]*")
        mult <- as.integer(str_extract(comb, "[[:digit:]].*"))
      }

      sum(round(PeriodicTable::mass(atom)) * mult)
    }))
  })
}
