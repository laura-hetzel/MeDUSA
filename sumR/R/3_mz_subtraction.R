# *** Subtraction -----------------------------------------------------
#' MZ-OBJ Subtraction (usually used for blanks subtraction)
#'
#' Compares the intensities of two "types" of samples.\cr
#'     Zeros any sample-intensities
#'     that are lower  than "Threshold"*"ComparedIntensities"
#'
#' @param sample_mz \cr
#'   DataFrame : Input MZ-Obj of Samples
#' @param subtract_mz \cr
#'   DataFrame : Metadata-Obj of Samples to subtract from
#' @param method: \cr
#'   String : method to apply to compare samples before threshold
#' @param threshold \cr
#'   Int       : "Delete samples < compare_samples * threshold"
#'
#' @return Returns an MZ-OBJ
#' @examples
#'
#' To remove values from "samples" that are lower than "blanks" * "threshold"
#'
#' sample_mz   <- sumR::mztools_filter(input_mzObj,metadata,"blanks","type",F)
#' subtract_mz <- sumR::mztools_filter(input_mzObj,metadata,"blanks","type",T)
#'
#' mz_subtraction(input_mz_obj, samples_meta, blanks_meta, threshold = "5" )
#
#' @export
mz_subtraction <- function(sample_mz_obj, subtract_mz_obj , method = mean, threshold = 3) {
  df_l <- local.ensure_mz(sample_mz_obj, subtract_mz_obj, "sumR::mz_subtraction")
  df_l$df_b["threshold"] <- apply(df_l$df_b, 1, method) * threshold

  .applyer <- function(input_mz_obj, compare = df_l$df_b["threshold"]){
    (input_mz_obj > compare) * input_mz_obj
  }
  df_sample <- as.data.frame(sapply(df_l$df_a, .applyer))
  df_sample <- cbind(df_l$mz, df_sample)
  out_mz <- df_sample[rowSums(dplyr::select(df_sample,-mz)) > 0,]
  local.mz_log_removed_rows(df_l$df_a,out_mz,"sumR::mz_subtraction")
  out_mz
}


