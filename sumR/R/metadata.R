#' @title Add sample data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with samplenames as rows
#' @export
addCellData <- function(exp, data = NULL){
  if (!validateExperiment(exp, checkColData = F)) return(NULL)

  if (is.null(data)) return(exp)
  cols <- rownames(colData(exp))
  data <- data[cols, , drop = F]

  columns <- c(colnames(colData(exp)), colnames(data))
  colData(exp) <- DataFrame(colData(exp), data)
  colnames(colData(exp)) <- columns
  exp
}

#' @title Add feature data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with featurenames as rows
#' @export
addFeatureData <- function(exp, data = NULL){
  if (!validateExperiment(exp)) return(NULL)

  if (is.null(data)) return(exp)
  data <- data[rownames(exp), ]
  rowData(exp) <- data.frame(rowData(exp), data)
  exp
}

#' @title Set metadata
#' @export
setMetadata <- function(exp, ...){
  if (!validateExperiment(exp)) return(NULL)

  metadata(exp) <- c(metadata(exp), list(...))
  exp
}

setDefaultAssay <- function(exp, default){
  if (!validateExperiment(exp)) return(NULL)
  n <- assayNames(exp)
  to_replace <- which(n == default)
  n[to_replace] <- n[1]
  n[1] <- default
  assays(exp) <- assays(exp)[n]
  exp
}

#' @title Set the phenotype column of a sumR Experiment
#' @param exp SummarizedExperiment obtained after alignment
#' @param phenotype Character vector specifying the column in colData to
#' be used as phenotype
#' @export
setPhenotype <- function(exp, phenotype){

  if (validateExperiment(exp)) {
    if (phenotype %in% colnames(colData(exp)) & is(phenotype, "character")) {
      metadata(exp)$phenotype <- phenotype
    } else {
      warning(sprintf("%s not found in colData", phenotype))
    }
  }
  exp
}
