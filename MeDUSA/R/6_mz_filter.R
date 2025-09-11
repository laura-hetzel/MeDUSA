# *** Blacklist -----------------------------------------------------
#' MZ-OBJ : Remove Blacklisted MZs
#'
#' @description
#' Filtering of the data is necessary to remove noise and contaminants from the
#' data set before statistical analysis. The mz_filter_blacklist function
#' accepts a list of m/z values that have been previously identified as
#' contaminants and removes them from the data set. This function can be used
#' to exclude any m/z values for any reason, but was designed for contaminants.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param blacklist \cr
#'   List/Vector : List of Blacklisted mzs, i.e. c(50.123,55.321)
#'   String      : Filename of Blacklisted mzs (csv: name,mass)
#'   Default     : Uses MeDUSA default_inputs
#' @param tolerance \cr
#'   Float   : Remove MZs +/- tolerance in ppm
#'
#' Dependencies :dplyr
#' @returns MZ-OBJ
#' @examples
#'   mz_filter_blacklist(input_mz, blacklist = c(100.111, 200.222), tolerance = 1e-7)
#'     : Remove blacklisted mzs within a given tolerance
#' @export
mz_filter_blacklist <- function( input_mz_obj, blacklist = "", tolerance = 5e-6){
  if ( class(blacklist) == "character"){
    if (blacklist == ""){
      blacklist <- get_default_data('blacklist')$value
    } else {
      blacklist <- read.csv(blacklist)$mass
    }
  }
  if (sum( blacklist < 0 ) > 0){
    stop("ERROR:MeDUSA::mz_filter_blacklist. Blacklist cannot have negative values.")
  } else if( length(blacklist) < 1 ){
    stop("ERROR:MeDUSA::mz_filter_blacklist. Blacklist is empty.")
  }
  bool_list <- lapply(blacklist, function(x){abs(input_mz_obj$mz - x)/x > tolerance})
  out_mz <- input_mz_obj[as.logical(Reduce("*",bool_list)),]
  local.mz_log_removed_rows(input_mz_obj, out_mz, "MeDUSA::mz_filter_blacklist")
  out_mz
}

# *** Misingness -----------------------------------------------------
#' MZ-OBJ or MZT-OBJ Filter Missingness
#
#' @description
#' Filtering of the data is necessary to remove noise and contaminants from the
#' data set before statistical analysis. The mz_filter_missingness function
#' removes m/z values that have a zero value in a specified percentage of the
#' samples.
#'s
#' Remove if peaks are significantly present.\cr
#'     Zeros any intensities that are not in [Threshold]% of samples
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param threshold \cr
#'   Numeric   : If Decimal : Threshold % to blank data under. (as a decimal 0.1 = 10%)
#'             : If Integer : Number of scans required to be nonzero ( 2 = Value required in at least 2 samples )
#' @param msg \cr
#'   String    : Identifier for log and plot outputs
#' @return Returns an MZ-OBJ
#' @examples
#'   mz_filter_missingness(input_mz_obj) : Run missingness filtering with custom values
#'   mz_filter_missingness(input_mz_obj, threshold = 3, msg = '3 thres')
#'     : Filter with custom missingness and logging message
#' @export
mz_filter_missingness <- function(input_mz_obj, threshold = 0.1, msg = ""){
  if ( threshold < 1){
    threshold <- threshold * (length(input_mz_obj)-1)
  }
  ##might be an issue with this line?
  input_mz_obj[is.na(input_mz_obj)] <- 0

  keep_peaks <- input_mz_obj[rowSums( dplyr::select(input_mz_obj, -mz) > 0 ) >= threshold,]
  local.mz_log_removed_rows(input_mz_obj,keep_peaks,paste("MeDUSA::mz_filter_missingness",msg))
  keep_peaks
}

# *** LowIntensity -----------------------------------------------------
#' MZ-OBJ or MZT-OBJ Filter LowIntensity
#'
#' @description
#' Filtering of the data is necessary to remove noise and contaminants from the
#' data set before statistical analysis. The mz_filter_lowIntensity function
#' will filter out noise by eliminating low intensity values. The user can set a
#' threshold for the minimum intensity, and any intensity below this value will
#' be set to zero.
#'
#' Remove if intensity is too low\cr
#'     Zeros any intensities that below [intensity].
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param threshold \cr
#'   Float     : Threshold to blank data under. (as an intensity.)
#'             : Suggested values: (Pos=10000, Neg=5000)
#' @param log_name \cr
#'   String    : Identifier for log and plot outputs
#' @return Returns an MZ-OBJ
#' @examples
#'   mz_filter_low_intensity(input_mz_obj) : Run low intensity filtering with custom values
#'   mz_filter_low_intensity(input_mz_obj, threshold = 5000, msg = '5000 thres')
#'     : Filter with custom intensity and logging message
#' @export
mz_filter_low_intensity <- function(input_mz_obj, threshold = 500, msg = ""){
  #TODO make this better
  mz <- input_mz_obj$mz
  input_mz_obj[is.na(input_mz_obj)] <- 0
  input_mz_obj[ input_mz_obj <= threshold ] <- 0
  input_mz_obj$mz <- mz
  out_mz <- input_mz_obj[ rowSums(dplyr::select(input_mz_obj ,-mz)) > 0, ]
  local.mz_log_removed_rows(input_mz_obj, out_mz, paste("MeDUSA::mz_filter_lowIntensity",msg))
  out_mz
}

# *** Filter-Magic -----------------------------------------------------
#' MZ-OBJ Filter-Magic
#'
#' @description
#' Filtering of the data is necessary to remove noise and contaminants from the
#' data set before statistical analysis. The mz_filter_magic function will
#' filter the data based on three parameters: the minimum intensity, a
#' missingness threshold, and a blacklist. The minimum intensity will set any
#' intensity value that is below the minimum to zero in an attempt to decrease
#' the statistical impact of noise. The missingness threshold will eliminate
#' m/z values that do not occur in a sufficient number of the samples. The
#' blacklist will eliminate m/z values that have been previously identified as
#' contaminants.
#'
#' Apply all filters as suggested by MeDUSA.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param min_intensity \cr
#'   Int       : Threshold to blank data under. (as a intesity value)\cr
#'             : Default values: (Pos=10000, Neg=5000)
#' @param missingness_thresholds \cr
#'   Float     : Threshold to blank data under. (as a decimal 10% = 0.1)
#' @param blacklist \cr
#'   List/Vector : List of Blacklisted mzs, i.e. c(50.123,55.321)
#'   String      : Filename of Blacklisted mzs (csv: name,mass)
#'   Default     : Uses MeDUSA default_inputs
#' Dependencies : dplyr
#' @return Returns an MZ-OBJ
#' @examples
#'   mz_filter_magic(input_mz_obj) : Run all filtering with default values
#'   mz_filter_magic(input_mz_obj, min_intensity = 5000, missingness_threshold = 2, blacklist = c(100.111, 200.222))
#'     : Run all filtering with custom values
#' @export
mz_filter_magic <- function(input_mz_obj, min_intensity = NULL, missingness_threshold = NULL, blacklist = NULL){
  if(!is.null(blacklist)){
    if( blacklist == TRUE ){
      input_mz_obj <- mz_filter_blacklist(input_mz_obj)
    } else {
      input_mz_obj <- mz_filter_blacklist(input_mz_obj, blacklist=blacklist)
    } 
  }
  if(!is.null(missingness_threshold)){
    input_mz_obj <- mz_filter_missingness(input_mz_obj, missingness_threshold)
  } else {
    default_missingness <- 0.1
    print(paste0("INFO:MeDUSA::mz_filter_magic: Using 'magic' missingness_threshold:",  default_missingness))
    input_mz_obj <- mz_filter_missingness(input_mz_obj, threshold = default_missingness )
  }
  if(is.null(min_intensity)){
    min_intensity <- tryCatch({
      local.mz_polarity_guesser(input_mz_obj, pos_return=5000, neg_return=2000)
    }, error = function(e) {
      print(e)
      stop("ERROR:MeDUSA::mz_filter_magic: Could not guess positive or negative from colnames
        Please specify min_intensity")
    })
    print(paste0("INFO:MeDUSA::mz_filter_magic: Using 'magic' min_intensity: ", min_intensity))
  }
  input_mz_obj <- mz_filter_low_intensity(input_mz_obj,min_intensity)
  input_mz_obj
}
