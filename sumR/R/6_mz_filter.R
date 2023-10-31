# *** Blacklist -----------------------------------------------------
#' MZ-OBJ Blacklist
#'
#' Removed Blacklisted MZs .\cr
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param blacklist \cr
#'   List/Vector : List of Blacklisted mzs, i.e. c(50.123,55.321)
#'   String      : Filename of Blacklisted mzs (csv: name,mass)
#'   Default     : Uses SumR default_inputs
#' @param tolerance \cr
#'   Float   : Remove MZs +/- tolerance in ppm
#'
#' Dependencies :dplyr
#' @return Returns an MZ-OBJ
#' @export
mz_filter_blacklist <- function( input_mz_obj, tolerance = 5e-6,
                        blacklist = sprintf("%s/../default_inputs/mz_blacklist.csv",getSrcDirectory(mz_tagging_blacklist))){
  if ( class(blacklist) == "character"){
    blacklist <- read.csv(blacklist)$mass
  }
  bool_list <- lapply(blacklist, function(x){abs(input_mz_obj$mz - x)/x > tolerance})
  out_mz <- input_mz_obj[as.logical(Reduce("*",bool_list)),]
  local.mz_log_removed_rows(input_mz_obj, out_mz, "sumR::mz_tagging_blacklist")
  out_mz
}

# *** Misingness -----------------------------------------------------
#' MZ-OBJ or MZT-OBJ Filter Missingness
#'
#' Remove if peaks are significantly present.\cr
#'     Zeros any intensities that are not in [Threshold]% of samples
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param threshold \cr
#'   Numeric   : If Decimal : Threshold % to blank data under. (as a decimal 10% = 0.1)
#'             : If Integer : Number of scans required to be nonzero ( i.e. Value required in at least 2 samples )
#' @return Returns an MZ-OBJ
#' @export
mz_filter_missingness <- function(input_mz_obj, threshold = 0.1, msg = ""){
  if ( threshold < 1){
    threshold <- threshold * (length(input_mz_obj)-1)
  }
  keep_peaks <- input_mz_obj[rowSums( dplyr::select(input_mz_obj, -mz) > 0 ) >= threshold,]
  local.mz_log_removed_rows(input_mz_obj,keep_peaks,paste("sumR::mz_filter_missingness",msg))
  keep_peaks
}

# *** LowIntensity -----------------------------------------------------
#' MZ-OBJ or MZT-OBJ Filter LowIntensity
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
mz_filter_lowIntensity <- function(input_mz_obj, threshold, msg = ""){
  #TODO make this better
  mz <- input_mz_obj$mz
  input_mz_obj[is.na(input_mz_obj)] <- 0
  input_mz_obj[ input_mz_obj <= threshold ] <- 0
  input_mz_obj$mz <- mz
  out_mz <- input_mz_obj[ rowSums(dplyr::select(input_mz_obj ,-mz)) > 0, ]
  local.mz_log_removed_rows(input_mz_obj, out_mz, paste("sumR::mz_filter_lowIntensity",msg))
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
#' @param missingness_thresholds \cr
#'   Float     : Threshold to blank data under. (as a decimal 10% = 0.1)
#'
#' Dependencies : dplyr
#' @return Returns an MZ-OBJ
#' @export
mz_filter_magic <- function(input_mz_obj, min_intensity, missingness_threshold=F, blacklist=T){
  if( sum(blacklist)){
    if (sum(blacklist) != T){
      input_mz_obj <- mz_filter_blacklist(input_mz_obj, blacklist = blacklist)
    }else{
      input_mz_obj <- mz_filter_blacklist(input_mz_obj)
    }
  }
  if(hasArg(missingness_threshold)){
    input_mz_obj <- mz_filter_missingness(input_mz_obj, missingness_threshold)
  } else {
    input_mz_obj <- mz_filter_missingness(input_mz_obj)
  }
  if(!hasArg(min_intensity)){
    min_intensity <- tryCatch({
      local.mz_polarity_guesser(input_mz_obj, pos_return=10000, neg_return=5000)
    }, error = function(e) {
      print(e)
      stop("ERROR: sumR::mz_filter_magic: Could not guess positive or negative from colnames
        Please specify min_intensity")
    })
  }
  input_mz_obj <- mz_filter_lowIntensity(input_mz_obj,min_intensity)
  input_mz_obj
}
