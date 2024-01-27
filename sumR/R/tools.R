
# *** Filter MZs by meta Phenotypes  -----------------------------------------------------
#' MZLOG-OBJ Plot Heat map
#'
#' Create Welch, Fold & Volcano plot, with default inputs
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param metadata \cr
#'   DataFrame : Metadata-Obj of Samples
#' @param phenotype \cr
#'   String: Phenotype to filter on
#' @param exclude \cr
#'.  Boolean: F = return values that do not match filter
#' @param keep_mz \cr
#'.  Boolean: should returned mzObj include mz columnd?
#'
#' @export
mztools_filter <- function(input_mzobj, metadata, filter_value , filter_name = "phenotype", exclude = F, keep_mz = T){
  if (exclude) {
    meta_filtered <- dplyr::filter(metadata, get(filter_name) != filter_value & filtered_out == "no")
  } else {
    meta_filtered <- dplyr::filter(metadata, get(filter_name) == filter_value & filtered_out == "no")
  }
  out <- input_mzobj[,meta_filtered$sample_name]
  if (keep_mz){
    out$mz <- input_mzobj$mz
    out <- dplyr::select(out,mz,everything())
  }
  
  out
}


