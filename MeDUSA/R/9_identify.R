# *** HMDB Identify -----------------------------------------------------
#' mz: Identify mzs in HMDB
#'
#' @description
#' A list of m/z values is not necessary helpful in metabolomic studies; these
#' values must be assigned to something meaningful. The identify_hmdb function
#' compares the m/z values previously identified as "relevant to phenotype
#' prediction" to the Human Metabolome Database (HMDB). Adducts may be added as
#' a list to improve the results of the search. Additionally, a tolerance is
#' set by the user for identifying matches with the database.
#'
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
identify_hmdb <- function( mzs, adducts = c("M+H"), hmdb_file = "/home/rstudio/local/hmdb_xml.Rdata", tolerance = 5e-6) {
  adducts <- identify.adducts(adducts, "identify_hmdb")
  #xmlToDataFrame takes ~15mim.
  if (grepl("Rdata$",hmdb_file)){
    load(hmdb_file)
  } else {
    hmdb <- XML::xmlToDataFrame(hmdb_file)
    #TODO save file?
    hmdb$monisotopic_molecular_weight <- as.numeric(hmdb$monisotopic_molecular_weight)
    hmdb <- hmdb[!is.na(hmdb$monisotopic_molecular_weight),]
  }
  mz_expanded <- unlist(lapply(mzs ,function(x) x - adducts))
  #TODO PBL apply
  out <- lapply(mz_expanded, function(x)
    hmdb[abs(hmdb$monisotopic_molecular_weight - x)/x < tolerance ,]
  )
  do.call(rbind,out)
}


# *** Lipid Identify -----------------------------------------------------
#' mz: Identify mzs in Lipids
#'
#' @description
#' A list of m/z values is not necessary helpful in lipidomic studies; these
#' values must be assigned to something meaningful. The identify_lipids function
#' compares the m/z values previously identified as "relevant to phenotype
#' prediction" to a lipids database. Adducts may be added as
#' a list to improve the results of the search. Additionally, a tolerance is
#' set by the user for identifying matches with the database.
#'
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
identify_lipids <- function( mzs, adducts = c("M+H"), lipids_file = "/home/rstudio/lipids_simple.csv", tolerance = 5e-6) {
  identify_from_csv(mzs, adducts, lipids_file, mz_colname = "EXACT_MASS", tolerance)
}

# *** Identify from csv-----------------------------------------------------
#' mz: Identify mzs in a provided csv
#'
#' @description
#' A list of m/z values is not necessary helpful in metabolomic studies; these
#' values must be assigned to something meaningful. It is possible to download
#' the desired databases for comparison and compare the data sets to the stored
#' CSV file with the identify_from_csv function. Adducts may be added as
#' a list to improve the results of the search. Additionally, a tolerance is
#' set by the user for identifying matches with the database.
#'
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
identify_from_csv <- function( mzs, adducts = c("M+H") , csv_file, mz_colname, tolerance = 5e-6) {
  adducts <- identify.adducts(adducts, "identify_lipids")

  if (grepl("Rdata$",csv_file)){
    load(csv_file)
  } else {
    csv_file <- readr::read_csv(csv_file)
    csv_file[[mz_colname]] <- as.numeric(csv_file[[mz_colname]])
    csv_file <- csv_file[!is.na(csv_file[[mz_colname]]),]
  }
  mz_expanded <- unlist(lapply(mzs ,function(x) x - adducts))
  out <- lapply(mz_expanded, function(x)
    csv_file[abs(csv_file[[mz_colname]] - x)/x < tolerance ,]
  )
  do.call(rbind,out)
}

identify.adducts <- function(adducts,source){
  if( is.vector(adducts)){
    if ( class(adducts) == "character"){
      common_adducts <- get_default_data('adducts')
      adducts = common_adducts$value[common_adducts$name %in% adducts]
    } else if ( class(adducts) != "numeric"){
      stop(paste0("ERROR: MeDUSA::identify::",source,"adducts cannot be:", class(adducts),". Examples of valid adducts: c(1.0008,17.03052) or c(\"H\")"))
    }
  } else {
    stop(paste0("ERROR: MeDUSA::identify::",source,"adducts is not a vector. Examples of valid adducts: c(1.0008,17.03052) or c(\"H\")"))
  }
  adducts
}
