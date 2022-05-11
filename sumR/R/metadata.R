#' @title Add sample data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with samplenames as rows
#' @export
addSampleData <- function(exp, data = NULL){
  if (is.null(data)) return(exp)
  cols <- rownames(colData(exp))
  data <- data[cols, ]
  colData(exp) <- cbind(colData(exp), data)
  exp
}

#' @title Add feature data to the experiment
#' @param exp SummarizedExperiment object
#' @param data Dataframe with featurenames as rows
#' @export
addFeatureData <- function(exp, data = NULL){
  if (is.null(data)) return(exp)
  data <- data[rownames(exp), ]
  rowData(exp) <- cbind(rowData(exp), data)
  exp
}


