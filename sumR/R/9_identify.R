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
identify_hmdb <- function( mzs, adduct, hmdb_file = "/home/rstudio/hmbd_metabolites.xml" ) {
  adducts <- c( "H"   = 1.0008 ,
                "NH4" = 17.03052  )
  hmdb <- XML::xmlParse("/home/rstudio/hmbd_metabolites.xml")

  #TODO, draw the rest of the owl

}
