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

#' @title Obtain metadata from aliquot names
#' @param files
#' @param regex
#' @export
metadataFromFile <- function(files, regex = "([0-9]{4}[A-Z]{3}_[0-9]{4})([A-Z]+|)([0-9]+|)_([A-Z])([0-9])"){
  aliquots <- tools::file_path_sans_ext(basename(files))
  matches <- regmatches(aliquots, gregexec(regex, aliquots))
  df <- do.call(rbind, lapply(matches, t))
  df <- as.data.frame(df)
  if (ncol(df) != 6) stop("Please check you regex or files, cannot parse")
  colnames(df) <- c("Aliquot", "Sample", "Type", "CalNo", "Replicate", "Injection")
  df$Type[df$Type == ""] <- "SAMPLE"
  df$CalNo[df$CalNo == ""] <- NA
  df$Replicate <- vapply(df$Replicate, utf8ToInt, numeric(1)) - 64
  df
}



