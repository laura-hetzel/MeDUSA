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
#'   String: Location of hmdb.xml, or hmdb.Rdata download file:
#'           DockerLocation: /home/rstudio/local/hmdb_simple.xml
#'
#' @export
identify_hmdb <- function( mzs, adducts = c("H"), hmdb_file = "/home/rstudio/hmdb_xml.Rdata", tolerance = 5e-6) {
  adducts <- .identify_adducts(adducts, "identify_hmdb")
  #xmlToDataFrame takes ~15mim.
  if (grepl("Rdata$",hmdb_file)){
    load(hmdb_file)
  } else {
    hmdb <- XML::xmlToDataFrame(hmdb_file)
    #TODO save file?
    hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
    hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  }
  mz_expanded = unlist(lapply(mzs ,function(x) x - adducts))
  print(mz_expanded)
  cl <- local.export_thread_env(2,environment())
  out <- lapply(mz_expanded, 1, function(x) 
    hmdb[abs(hmdb$monisotopic_molecular_weight - x)/x < tolerance ,]
  )
  local.kill_threads(cl)
  do.call(rbind,out)
}


# *** Lipid Identify -----------------------------------------------------
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
#' @param lipids \cr
#'   String: Location of lipids.csv, or hmdb.Rdata download file:
#'           DockerLocation: /home/rstudio/local/lipids_simple.csv
#'
#' @export
identify_lipids <- function( mzs, adducts = c("H") , lipids_file = "/home/rstudio/lipids_simple.csv", tolerance = 5e-6) {
  adducts <- .identify_adducts(adducts, "identify_lipids")

  if (grepl("Rdata$",lipids_file)){
    load(lipids_file)
  } else {
    lipids <- readr::read_csv(lipids_file)
    lipids$EXACT_MASS <- as.numeric(lipids$EXACT_MASS)
    lipids <- lipids[!is.na(lipids$EXACT_MASS),]
  }
  mz_expanded = unlist(lapply(mzs ,function(x) x - adducts))
  print(mz_expanded)
  #cl <- local.export_thread_env(2,environment())
  out <- lapply(mz_expanded, function(x) 
    lipids[abs(lipids$EXACT_MASS - x)/x < tolerance ,]
  )
  #local.kill_threads(cl)
  do.call(rbind,out)
}


.identify_adducts <- function(adducts,source){
  if( is.vector(adducts)){
    if ( class(adducts) == "character"){
      common_adducts <- data.frame( "name" = c("H","NH4","K"),
                                    "mz"   = c(1.0008, 17.03052 ,39.0983 ))
      adducts = common_adducts$mz[common_adducts$name %in% adducts]
    } else if ( class(adducts) != "numeric"){
      stop(paste0("ERROR: sumR::identify::",source,"adducts cannot be:", class(adducts),". Examples of valid adducts: c(1.0008,17.03052) or c(\"H\")"))
    }
  } else {
    stop(paste0("ERROR: sumR::identify::",source,"adducts is not a vector. Examples of valid adducts: c(1.0008,17.03052) or c(\"H\")"))
  } 
  adducts
}