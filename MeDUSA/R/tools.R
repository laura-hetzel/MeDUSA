# *** Filter MZs by meta Phenotypes  -----------------------------------------------------
#' MzObj: filter out samples via meta queries
#'
#' @description
#' The mztools_filter function is used to filter the m/z values by phenotype.
#' This filtering is useful for t-test, fold change, and volcano plots. For
#' optimal results, the metadata should have a properly populated "phenotype"
#' column.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param metadata \cr
#'   DataFrame : Metadata-Obj of Samples
#' @param filter_value \cr
#'   String: Value of the metadata column to filter
#' @param filter_name \cr
#'   String: Name of the metadata column to filter
#' @param exclude \cr
#'.  Boolean: F = return values that do not match filter
#' @param keep_mz \cr
#'.  Boolean: should returned mzObj include mz columnd?
#'
#' @export
mztools_filter <- function(input_mzobj, metadata, filter_value , filter_name = "phenotype", exclude = F, keep_mz = T){
  metadata <- local.meta_polarity_fixer(input_mzobj, metadata)
  if (exclude) {
    meta_filtered <- dplyr::filter(metadata, !(get(filter_name) %in% filter_value) & filtered_out == FALSE)
  } else {
    meta_filtered <- dplyr::filter(metadata, get(filter_name) %in% filter_value & filtered_out == FALSE)
  }
  out <- input_mzobj[,meta_filtered$sample_name]
  if (keep_mz){
    out$mz <- input_mzobj$mz
    out <- dplyr::select(out,mz,everything())
  }
  out
}

# *** Get Default Data-----------------------------------------------------
#'
#' The get_default_data function will return a list of default values including
#' available adducts and their corresponding masses as well as default blacklist
#' of features to exclude and their corresponding masses.
#'
#' @param type \cr
#'   Character : 'adduct' or 'blacklist'
#'
#' @returns dataframe( [adduct]name: Character,
#'                     [adduct]value: Numeric )
#'
#' @export
get_default_data <- function(type){
  if( type == 'adducts'){
    data_list <- list(
      c( "M+H",   -1.0008   ),
      c( "M-H",   +1.0008   ),
      c( "M+NH4", -14.0067  ),
      c( "2M+H",  -0.5004   ),
      c( "2M-H",  +0.5004   ),
      c( "M+Na",  -22.990   ),
      c( "M+Cu",  -63.546   ),
      c( "2M+Na", -11.4950  ),
      c( "2M+Cu", -31.773   ),
      c( "3M+Na", -7.663333 ),
      c( "CN-",   -16.018   ),
      c( "HCOO-", -46.0254  ),
      c( "M+FormicAcid+H", -47.0262 ),
      c( "M+DimethylFormamide+H", -74.0946 )
    )
  } else if (type == 'blacklist') {
    data_list <- list(
      c("Polysiloxane", 536.17)
    )
  } else {
    stop("ERROR: MeDUSA::get_default_values: Data type not found.")
  }
  out <- data.frame(name = character(), value = numeric())
  for(x in 1:length(data_list)){
    out <- tibble::add_row(out, name=data_list[[x]][1], value=as.numeric(data_list[[x]][2]))
  }
  out
}
