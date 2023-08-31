# *** Subtraction -----------------------------------------------------
#' MZ-OBJ Subtraction (usually used for blanks subtraction)
#'
#' Compares the intensities of two "types" of samples.\cr
#'     Zeros any sample-intensities
#'     that are lower  than "Threshold"*"ComparedIntensities"
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param samples \cr
#'   DataFrame : Metadata-Obj of Samples
#' @param compare_samples \cr
#'   DataFrame : Metadata-Obj of Compared samples (blanks)
#' @param threshold \cr
#'   Int       : "Delete samples < compare_samples * threshold"
#'
#' @return Returns an MZ-OBJ
#' @examples
#'
#' To remove values from "samples" that are lower than "blanks" * "threshold"
#'
#'
#' blanks_meta <- filter(metadata, sampletype == "blank" &
#'                                   ionmode == ionmode_val & phase == phase_val &
#'                                   filtered_out == "no")
#' samples_meta <- filter(metadata, sampletype == "sample" &
#'                                   ionmode == ionmode_val  & phase == phase_val &
#'                                   filtered_out == "no")
#'
#' mz_subtraction(input_mz_obj, samples_meta, blanks_meta, threshold = "5" )
#
#' @export

mz_subtraction <- function(input_mz_obj, samples, compare_samples , threshold = 3) {
  input_mz_obj <- as.data.frame(input_mz_obj)

  compare_subset <- input_mz_obj[,compare_samples$sample_name]
  compare_subset[threshold] <- apply(compare_subset, 1, median) * threshold

  samples_subset <- input_mz_obj[,samples$sample_name]

  .applyer <- function(input_mz_obj, compare = compare_subset[threshold]){
    (input_mz_obj > compare) * input_mz_obj
  }
  samples_subset <- as.data.frame(sapply(samples_subset, .applyer))
  row.names(samples_subset) <- row.names(input_mz_obj)
  out_mz <- samples_subset[rowSums(samples_subset) > 0,]
  local.mz_log_removed_rows(input_mz_obj,out_mz,"sqrlSumr::mz_subtraction")
  out_mz
}
