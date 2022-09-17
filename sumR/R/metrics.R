#' @title Add the Sample/Blank ratio for each compound
#' @param exp SummarizedExperiment Object
#' @param typeColumn Name of the column in `colData(exp)` that represents
#' the sample type and contains at least the values `SAMPLE` and `BLANK`.
#' @export
addSampleBlankRatio <- function(exp, assay = 1, typeColumn = "Type"){
  samps <- as.matrix(assay(exp[, exp$Type == "SAMPLE"], assay))
  blanks <- as.matrix(assay(exp[, exp$Type == "BLANK"], assay))
  rowData(exp)$SampleBlankRatio <- round(rowMedians(samps, na.rm = TRUE) /
                                           rowMedians(blanks, na.rm = TRUE), 3)
  exp
}
