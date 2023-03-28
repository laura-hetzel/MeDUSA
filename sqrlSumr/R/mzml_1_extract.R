# ***  -----------------------------------------------------
#' MZML Extract
#'
#' Not really sure
#'  - Requires: Bioc::mzR
#'
#' @param input_mzml \cr
#'   String : Filename
#' @returns [MZ-OBJ,MZ-Meta]
#'
#' @export
mzml_extract <- function(input_mzml) {
  mzML_xml <- mzR::openMSfile(input_mzml,  backend = "pwiz")
  mzML_header <- mzR::header(mzML_xml)
  isPos <- grepl("\\+",mzML_header$filterString)]
  mzML_negScans <- mzML_header[!isPos,]$seqNum
  mzML_posScans <- mzML_header[isPos,]$seqNum
  ]
  mzML_spectra_pos <- mzR::spectra(mzML_xml,mzML_posScans)
}
