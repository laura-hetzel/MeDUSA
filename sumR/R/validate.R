validateSpectra <- function(spectraDf){
  tests <- c(
    class(spectraDf) == "DFrame",
    attr(spectraDf, "package") == "S4Vectors",
    all(c("sample", "mz", "mzmin", "mzmax", "npeaks", "i", "fpeaks",
          "ppm", "peakIdx", "rt") %in% colnames(spectraDf)),
    typeof(spectraDf$sample) == "character",
    typeof(spectraDf$mz) == "double",
    typeof(spectraDf$mzmin) == "double",
    typeof(spectraDf$mzmax) == "double",
    typeof(spectraDf$npeaks) == "integer",
    typeof(spectraDf$i) == "double",
    typeof(spectraDf$fpeaks) == "double",
    typeof(spectraDf$ppm) == "double",
    typeof(spectraDf$peakIdx) == "list",
    typeof(spectraDf$peakIdx[[1]]) == "character",
    typeof(spectraDf$rt) == "double",
    nrow(spectraDf) > 0
  )
  all(tests)
}


#' @title Validate if a centroided data.frame is valid to use with binning
#' @param peakDf data.frame with scan, mz, i, and Noise columns.
#' @noRd
validatePeak <- function(peakDf){
  att <- attributes(peakDf)
  all(
    class(peakDf) == "data.frame",
    all(c("scan", "mz", "i") %in% colnames(peakDf)),
    nrow(peakDf) > 0,
    all(c("polarity", "massWindow", "combineSpectra", "Files") %in% names(att)),
    any(c("-", "+") %in% att$polarity),
    any(c(TRUE, FALSE) %in% att$combineSpectra),
    length(att$Files) == 1
  )
}

#' @title Validate if a list of centroided data.frames is valid to
#' use with binning
#' @param peakList list of data.frames, each with scan, mz, i, and
#' Noise columns.
#' @noRd
validatePeaks <- function(peakList){
  all(
    length(names(peakList)) == length(peakList),
    all(sapply(peakList, validatePeak))
  )
}


#' @title Validate a sumR Experiment on rows, columns and assays
#' @param exp SummarizedExperiment obtained after alignment
#' @importFrom SummarizedExperiment assayNames
#' @export
validateExperiment <- function(exp, checkColData = T){
  isValid <- all(
    class(exp) == "SummarizedExperiment",
    nrow(exp) > 0,
    ncol(exp) > 0,
    #ifelse(checkColData, "Type" %in% colnames(colData(exp)), TRUE),
    #c("mz", "mzmin", "mzmax", "npeaks", "peakidx") %in% colnames(rowData(exp)),
    "Area" %in% assayNames(exp)
  )
  if (!isValid) warning("Invalid sumR Experiment object")
  isValid
}
