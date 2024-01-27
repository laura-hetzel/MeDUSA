# *** HMDB Identify -----------------------------------------------------
#' MZLOG-OBJ PCA
#'
#' Not really sure
#'  - Requires: ggplot2, tibble
#'
#' @param mzs \cr
#'   List : List of MZs to Identify
#' @param adduct \cr
#'   String: Name of a common adduct
#'   Integer: Atomic mass of an adduct
#' @param hmdb_file \cr
#'   String: Location of hmdb.xml download file
#'
#' @export
identify_hmdb <- function( mzs, adducts = 0 , hmdb_file = "/home/rstudio/local/hmdb_xml.Rdata", tolerance = 5e-6) {
  if (class(adducts) == "character"){
    common_adducts <- data.frame( "name" = c("H","NH4","K"),
                                  "mz"   = c(1.0008, 17.03052 ,39.0983 ))
    adducts = common_adducts$mz[common_adducts$name %in% adducts]
  }
  #xmlToDataFrame takes ~15mim.
  if (grep("Rdata$",hmdb_file)){
    load(hmdb_file)
  } else {
    hmdb <- XML::xmlToDataFrame(hmdb_file)
    hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
    hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  }
  mz_expanded = unlist(lapply(mzs ,function(x) x - adducts))
  print(mz_expanded)
  cl <- local.export_thread_env(2,environment())
  out <- lapply(mz_expanded, function(x) 
    hmdb[abs(hmdb$monisotopic_molecular_weight - x)/x < tolerance ,]
  )
  local.kill_threads(cl)
  do.call(rbind,out)
}


# *** HMDB Identify -----------------------------------------------------
#' MZLOG-OBJ PCA
#'
#' Not really sure
#'  - Requires: ggplot2, tibble
#'
#' @param mzs \cr
#'   List : List of MZs to Identify
#' @param adduct \cr
#'   String: Name of a common adduct
#'   Integer: Atomic mass of an adduct
#' @param hmdb_file \cr
#'   String: Location of hmdb.xml download file
#'
#' @export
identify_lipids <- function( mzs, adducts = c(0) , lipids_file = "/home/rstudio/local/hmdb_xml.Rdata", tolerance = 5e-6) {

  #xmlToDataFrame takes ~15mim.
  if (grep("Rdata$",hmdb_file)){
    load(hmdb_file)
  } else {
    hmdb <- XML::xmlToDataFrame(hmdb_file)
    hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
    hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  }
  mz_expanded = unlist(lapply(mzs ,function(x) x - adducts))
  print(mz_expanded)
  cl <- local.export_thread_env(2,environment())
  out <- lapply(mz_expanded, function(x) 
    hmdb[abs(hmdb$monisotopic_molecular_weight - x)/x < tolerance ,]
  )
  local.kill_threads(cl)
  do.call(rbind,out)
}


.identify_adducts <- function(adducts,source){
  if( is.vector(adducts)){
    if ( class(adducts) == "character"){
      common_adducts <- data.frame( "name" = c("H","NH4","K"),
                                    "mz"   = c(1.0008, 17.03052 ,39.0983 ))
      adducts = common_adducts$mz[common_adducts$name %in% adducts]
    } 
  } else {
    stop(paste0("Error: ",source,"adducts is not a vector. Examples of valid adducts: c(1.0008,17.03052) or c(\"H\")"))
  }
  
  
}