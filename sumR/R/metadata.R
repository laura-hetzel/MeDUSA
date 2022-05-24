#' @title Add sample data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with samplenames as rows
#' @export
addCellData <- function(exp, data = NULL){
  if (is.null(data)) return(exp)
  cols <- rownames(colData(exp))
  data <- data[cols, ]
  colData(exp) <- data.frame(colData(exp), data)
  exp
}

#' @title Add feature data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with featurenames as rows
#' @export
addFeatureData <- function(exp, data = NULL){
  if (is.null(data)) return(exp)
  data <- data[rownames(exp), ]
  rowData(exp) <- data.frame(rowData(exp), data)
  exp
}

#' @title Set metadata
#' @export
setMetadata <- function(exp, ...){
  metadata(exp) <- c(metadata(exp), list(...))
  exp
}

#' @title Create metadata from Excel file
#' @param xlsxFile
#' @param idxColumn
#' @param sheet
#' @export
metadataFromExcel <- function(xlsxFile, idxColumn, sheet = 1){
  df <- openxlsx::read.xlsx(xlsxFile, sheet)
  if (!idxColumn %in% colnames(df)) stop(sprintf("Error reading excel file, could not find column %s", idxColumn))
  if (any(duplicated(df[, idxColumn]))) stop(sprintf("Non-unique values found in %s", idxColumn))

  rownames(df) <- tools::file_path_sans_ext(basename(df[, idxColumn]))
  df
}


