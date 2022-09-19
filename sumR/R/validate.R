#' @title Validate the spectra DataFrame
#' @description This function validates if the DataFrame from `spectraAlignment`
#' is suitable for further processing.
#' @details `spectraAlignment` returns a DataFrame from S4Vectors containing the
#' following columns:
#'
#' * __sample__: Name of the sample the peak originated from
#' * __mz__: Median mz of the found peak
#' * __mzmin__: Lowest mz of the peak
#' * __mzmax__: Highest mz of the peak
#' * __npeaks__: Number of scans the peak consists of
#' * __i__: (Integrated) Intensity value of the peak
#' * __fpeaks__: Fraction of scans the peak is found in. Calculated by dividing
#' _npeaks_ by the total number of scans in the mz region.
#' * __ppm__: Parts per million error of the peak. Calculated by
#' \deqn{ \frac{mzmax - mzmin}{mz} \cdot 1e^6 }.
#' * __peakIdx__: List of indexes where the scans are found.
#' * __rt__: Median rt of the peak
#' @returns Either TRUE or FALSE representing the validation of the object.
#' @param spectraDf DataFrame object obtained from `spectraAlignment`
#' @importFrom methods is
#' @export
validateSpectra <- function(spectraDf){
  tests <- c(
    is(spectraDf, "DFrame"),
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

#' @title Validate if a list of centroided data.frames is valid to
#' use with alignment
#' @description This function checks if the centroiding went correctly and
#' the list of data.frames contain the right columns to be used for
#' alignment. See details about the validation procedure.
#' @details This function checks the list of data.frames on the following
#' points:
#'
#' - Do the names of the list equal the number of entries in the list
#' - does each entry is a data.frames with columns `scan`, `rt`, `mz` and `i`.
#' - Does each entry have the attributes `polarity`, `massWindow`,
#' `combineSpectra` and `Files`.
#' @param peakList list of data.frames, each with scan, mz, i, and
#' rt columns.
#' @export
validatePeaks <- function(peakList){
  all(
    length(names(peakList)) == length(peakList),
    all(sapply(peakList, function(peakDf){
      att <- attributes(peakDf)
      all(
        class(peakDf) == "data.frame",
        all(c("scan", "mz", "i", "rt") %in% colnames(peakDf)),
        nrow(peakDf) > 0,
        all(c("polarity", "massWindow", "combineSpectra", "Files") %in% names(att)),
        any(c("-", "+") %in% att$polarity),
        any(c(TRUE, FALSE) %in% att$combineSpectra),
        length(att$Files) == 1
      )
    }))
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
    "Area" %in% assayNames(exp)
  )
  if (!isValid) warning("Invalid sumR Experiment object")
  isValid
}
