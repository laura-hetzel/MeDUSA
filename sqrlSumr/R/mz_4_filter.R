# *** Misingness -----------------------------------------------------
#' MZ-OBJ Filter Missingness
#'
#' Remove if peaks are significantly present.\cr
#'     Zeros any intensities that are not in [Threshold]% of samples
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param threshold \cr
#'   Float     : Threshold to blank data under. (as a decimal 10% = 0.1)
#'
#' @return Returns an MZ-OBJ
#' @export
mz_filter_missingness <- function(input_mz_obj, threshold = 0.1){
  thresh <- threshold * length(input_mz_obj)
  keep_peaks <- input_mz_obj[rowSums( input_mz_obj > 0 ) > thresh,]
  local.mz_log_removed_rows(input_mz_obj,keep_peaks,"sqrlSumr::mz_filter_missingness")
  keep_peaks
}

# *** LowIntensity -----------------------------------------------------
#' MZ-OBJ Filter LowIntensity
#'
#' Remove if intensity is too low\cr
#'     Zeros any intensities that below [intensity].
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param threshold \cr
#'   Float     : Threshold to blank data under. (as an intensity.)
#'             : Suggested values: (Pos=10000, Neg=5000)
#'
#' @return Returns an MZ-OBJ
#' @export
mz_filter_lowIntensity <- function(input_mz_obj, threshold){
  out_mz <- input_mz_obj[rowSums(input_mz_obj > threshold)>0,]
  local.mz_log_removed_rows(input_mz_obj, out_mz, "sqrlSumr::mz_filter_lowIntensity")
  out_mz
}

# *** Filter-Magic -----------------------------------------------------
#' MZ-OBJ Filter-Magic
#'
#' Apply all filters as suggested by SumR.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param min_intensity \cr
#'   Int       : Threshold to blank data under. (as a intesity value)\cr
#'             : Default values: (Pos=10000, Neg=5000)
#' @param min_missingness \cr
#'   Float     : Threshold to blank data under. (as a decimal 10% = 0.1)
#'
#' Dependencies : dplyr
#' @return Returns an MZ-OBJ
#' @export
mz_filter_magic <- function(input_mz_obj, min_intensity, min_missingness=FALSE){
  if(hasArg(min_missingness)){
    tmp_mz <- mz_filter_missingness(input_mz_obj,min_missingness)
  } else {
    tmp_mz <- mz_filter_missingness(input_mz_obj)
  }
  if(!hasArg(min_intensity)){
    min_intensity <- tryCatch({
      local.mz_polarity_guesser(input_mz_obj, pos_return=10000, neg_return=5000)
    }, error = function(e) {
      print(e)
      stop("ERROR: sqrlSumr::mz_filter_magic: Could not guess positive or negative from colnames
        Please specify min_intensity")
    })
  }
  tmp_mz <- mz_filter_lowIntensity(tmp_mz,min_intensity)
  tmp_mz
}
