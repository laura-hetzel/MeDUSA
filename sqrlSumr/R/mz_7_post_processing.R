# *** Imputation -----------------------------------------------------
#' MZ-OBJ Imputation
#'
#' Replace all <[low_noise] ("mostly zeros") with a randomized value
#'   between [low_noise] and [high_noise]
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param low_noise \cr
#'   Float     : LowBoundry for "Noise"
#' @param high_noise \cr
#'   Float     : HighBoundry for "Noise"
#'
#' @returns MZ-OBJ
#' @export
mz_pp_imputation <- function(input_mz_obj, low_noise=10, high_noise=5000){
  if(!hasArg(high_noise)){
    high_noise <- tryCatch({
      local.mz_polarity_guesser(input_mz_obj, pos_return=10000, neg_return=5000)
    }, error = function(e) {
      print(e)
      stop("ERROR: sqrlSumr::mz_pp_imputation: Could not guess positive or negative from colnames
            Please specify min_intensity")
    })
  }
  if (low_noise > 100 ){
    warning("WARN: sqrlSumr::mz_PP_imputation: low_noise is < 100.
      Make sure [input_mz_obj] does not have an mz column ")
  }

  mz_bool <- input_mz_obj < low_noise
  input_mz_obj[mz_bool] <- sample(low_noise:high_noise,
                                              sum(mz_bool),
                                              replace = TRUE)
  input_mz_obj
}

# *** Normalization -----------------------------------------------------
#' MZ-OBJ Normalization
#'
#' Not really sure what this does
#'  - Requires: ggplot2, dplyr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#'
#' @returns Returns an MZ-OBJ
#' @export
mz_pp_normalization <- function(input_mz_obj, metadata, plot = TRUE ){
  #TODO revisit, for a refactor
  normal <- quotNorm(t(input_mz_obj))
  dilution <- data.frame(normal$dilution)
  dilution$sample <- rownames(dilution)
  join_by <- c('sample')
  names(join_by) <- "sample_name"
  dilution <- dplyr::left_join(metadata, dilution, by = join_by)

  if(plot){
    ggplot(dilution,
           aes(x = phenotype, y = normal.dilution)) + geom_boxplot()
    local.save_plot(paste("Normalization",local.mz_polarity_guesser(input_mz_obj),sep="-"))
  }

  normal <- data.frame(normal$X)
  normal <- as.data.frame(t(normal))
  rownames(normal) <- gsub("X", "",rownames(normal)) #TODO fix 'quotNorm' to not return mz = "X100.2"
  normal
}

# *** PivotLonger -----------------------------------------------------
#' MZ-OBJ Pivot Longer
#'
#' Not really sure what this does
#'  - Requires: ggplot2, tidyr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#'
#' @returns MZLONG-OBJ
#' @export
mz_pp_pivot_longer <- function(input_mz_obj, plot = TRUE) {
  input_mz_obj$mz <- row.names(input_mz_obj)
  input_mzlong <- input_mz_obj %>%
    tidyr::pivot_longer(!mz, names_to = "sample", values_to = "intensity")
    if(plot){
      ggplot(input_mzlong, aes(x = sample, y = intensity)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(paste("Normalized, Filtered,",
          local.mz_polarity_guesser(input_mz_obj)
          , sep=""))
      local.save_plot(paste("PivotLonger",local.mz_polarity_guesser(input_mz_obj),sep="-"))
    }
  data.frame(input_mzlong)
}

# *** Log Transform -----------------------------------------------------
#' MZLONG-OBJ Log Transform
#'
#' Convert intensities of MZLong to log2
#'  For normal "mz_obj" use log2(mz_obj)
#'  - Requires: ggplot2
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZLONG-Obj
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#'
#' @export
mzl_pp_log_transform <- function(input_mzlong, plot = TRUE){
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

# *** Process Magic -----------------------------------------------------
#' MZ-OBJ Process Magic
#'
#' Do all the recommended post processing
#'  - Requires: ggplot2, tidyrx
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param metadata \cr
#'   DataFrame : Input MZ-Metadata
#' @param plot \cr
#'   Boolean   : to plot or not to plot
#'
#' @returns c(mz_long, mz_obj)
#' @export
mz_post_process_magic <- function(input_mz_obj, metadata, plot = TRUE){
  tryCatch({
    input_mz_obj <- mz_pp_imputation(input_mz_obj)
    print("INFO: imputation success")
  }, error = function(e) {
    print(e)
    stop("ERROR: in mz_pp_imputation")
  })
  tryCatch({
    input_mz_obj <- mz_pp_normalization(input_mz_obj, metadata, plot)
    print("INFO: normalization success")
  }, error = function(e) {
    print(e)
    stop("ERROR: in mz_pp_normalization")
  })
  tryCatch({
    input_mzlong <- mz_pp_pivot_longer(input_mz_obj, plot)
    print("INFO: pivot_longer success")
  }, error = function(e) {
    print(e)
    stop("ERROR: in mz_pp_pivot_longer")
  })
  tryCatch({
    input_mzlong <- mzl_pp_log_transform(input_mzlong, plot)
    input_mz_obj <- log2(input_mz_obj)
    print("INFO: log_transform success")
  }, error = function(e) {
    print(e)
    stop("ERROR: in mzl_pp_log_transform")
  })
  list(input_mzlong ,input_mz_obj)
}
