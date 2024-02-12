# *** HMDB Identify -----------------------------------------------------
#' MZLOG-OBJ PCA
#'
#' Not really sure
#'  - Requires: ggplot2, tibble
#'
#' @param mzs \cr
#'   List : List of MZs to Identify
#' @param polarity \cr
#'   String: 'pos' | 'neg'
#' @param adduct \cr
#'   String: Name of a common adduct
#'   Integer: Atomic mass of an adduct
#' @param hmdb_file \cr
#'   String: Location of hmdb.xml, or hmdb.Rdata download file:
#'           DockerLocation: /home/rstudio/local/hmdb_simple.xml
#'
#' @export
identify_hmdb <- function( mzs, polarity, adducts = c("H"),hmdb_file = "/home/rstudio/hmdb_xml.Rdata", tolerance = 5e-6) {
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
  mz_expanded = .identify_adduct_math(mz_expanded, adducts, polarity)
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
identify_lipids <- function( mzs, adducts = c("H"), polarity, lipids_file = "/home/rstudio/lipids_simple.csv", tolerance = 5e-6) {
  identify_from_csv(mzs, adducts, csv_file = lipids_file, tolerance)
}

# *** Identify from csv-----------------------------------------------------
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
identify_from_csv <- function( mzs, polarity, adducts = c("H") , csv_file, tolerance = 5e-6) {
  adducts <- .identify_adducts(adducts, "identify_lipids")

  if (grepl("Rdata$",csv_file)){
    load(csv_file)
  } else {
    csv_file <- readr::read_csv(csv_file)
    csv_file$EXACT_MASS <- as.numeric(csv_file$EXACT_MASS)
    csv_file <- csv_file[!is.na(csv_file$EXACT_MASS),]
  }
  mz_expanded = .identify_adduct_math(mz_expanded, adducts, polarity)
  cl <- local.export_thread_env(2,environment())
  out <- lapply(mz_expanded, function(x) 
    csv_file[abs(csv_file$EXACT_MASS - x)/x < tolerance ,]
  )
  local.kill_threads(cl)
  do.call(rbind,out)
}

.identify_adduct_math <- function(mzs, adducts, polarity){
  if ( polarity == "neg"){
    mz_expanded = unlist(lapply(mzs ,function(x) x - adducts))
  } else if (polarity == "pos"){
    mz_expanded = unlist(lapply(mzs ,function(x) x + adducts))
  } else {
    stop("sumR:: identify: Polarity not found. Use 'pos' or 'neg'")
  }
  mz_expanded
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