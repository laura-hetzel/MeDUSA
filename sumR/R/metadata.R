#' @title Add cell/sample metadata to the SummarizedExperiment object
#' @description This function adds the given metadata in `data` to the
#' [colData] slot of the SummarizedExperiment `exp`.
#' @details Metadata about cells/samples is vital for post-processing. For
#' each cell, its type (SAMPLE, BLANK) must be specified in the column `Type`
#' in order to perform blank filtering. Another column specifiying the
#' phenotype of a cell must also be present for statistics and/or modelling. By
#' default, the `Datetime` of the sample is extracted from the mzML data using
#' the [extractPeaks] function. The cellData can also be visualized with
#' various plots implented in sumR, e.g. [pheatmapPlot].
#'
#' After setting the cellData, the results can be seen in the [colData] slot
#' of the SummarizedExperiment.
#' @returns SummarizedExperiment with an updated [colData] slot.
#' @param exp SummarizedExperiment object obtained after alignment.
#' @param data Dataframe with samplenames as rows
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Read the cellData
#' cellData <- read.csv(system.file("cellData.csv", package = "sumR"),
#'                      sep = ",", row.names = 1)
#'
#' # Add the cellData
#' sumRnegative <- addCellData(sumRnegative, cellData)
#'
#' # Show the colData
#' colData(sumRnegative)
addCellData <- function(exp, data = NULL){
  if (!validateExperiment(exp)) return(exp)

  if (is.null(data)) return(exp)
  cols <- rownames(colData(exp))
  data <- data[cols, , drop = F]

  columns <- c(colnames(colData(exp)), colnames(data))
  colData(exp) <- DataFrame(colData(exp), data)
  colnames(colData(exp)) <- columns
  exp
}

#' @title Add metadata of peaks to the SummarizedExperiment
#' @description Peaks are the rows in a SummarizedExperiment. This function
#' allows you to set metadata about these peaks. In most cases, this is handled
#' by sumR itself, however, any additional annotation data may be set using
#' this function.
#' @details Metadata of peaks are stored in the rowData of the
#' SummarizedExperiment. These may contain metrics of peaks like mz and rt or
#' other summarized information like RSD, or other ratios. Therefore, rowData
#' is often set by sumR itself. In some cases, like compound annotation,
#' identifiers or other metadata it may be set by the user. This function aids
#' in doing so.
#' @returns SummarizedExperiment object with added peak data in the [rowData]
#' slot of the SummarizedExperiment
#' @param exp SummarizedExperiment object obtained after alignment
#' @param data Dataframe with peaknames as rows and attributes as columns.
#' @importFrom SummarizedExperiment rowData rowData<-
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Setup data.frame of polarity
#' df <- data.frame(Polarity = rep("positive",
#'                  nrow(sumRnegative)),
#'                  row.names = rownames(sumRnegative))
#'
#' # Adding polarity to the peaks
#' sumRnegative <- addPeakData(sumRnegative, df)
#'
#' # Showing rowData
#' rowData(sumRnegative)
addPeakData <- function(exp, data = NULL){
  if (!validateExperiment(exp)) return(NULL)

  if (is.null(data)) return(exp)

  rows <- rownames(rowData(exp))
  data <- data[rows, , drop = F]

  columns <- c(colnames(rowData(exp)), colnames(data))
  rowData(exp) <- DataFrame(rowData(exp), data)
  colnames(rowData(exp)) <- columns
  exp
}

#' @title Add metadata to the SummarizedExperiment
#' @description the [metadata] slot of the SummarizedExperiment can hold any
#' type data and is often used for describing the experiment or data types that
#' don't fit in the [assay], [colData], or [rowData] slots. This function aids
#' in setting metadata.
#' @details Metadata in a SummarizedExperiment is used for miscelaneous data
#' that describe the experiment. Examples are the datetime of when the
#' experiment was made or model information (which is what it is used for in
#' sumR). However, additional data can be set as well with this function.
#' @returns SummarizedExperiment object with added metadata in the [metadata]
#' slot.
#' @param exp SummarizedExperiment object obtained after alignment
#' @param metaList Named list of metadata fields that can be of any type.
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Set the current time in the experiment
#' sumRnegative <- addMetadata(sumRnegative, list(Time = lubridate::now()))
#'
#' # Show the metadata slot
#' metadata(sumRnegative)
addMetadata <- function(exp, metaList){
  if (!validateExperiment(exp)) return(NULL)

  metadata(exp) <- c(metadata(exp), metaList)
  exp
}

#' @title Set the default assay in the SummarizedExperiment
#' @description This function swaps the current assay at index 1  in `exp`
#' for the assay provided in `default`.
#' @details By default, most post-processing functions use the assay in `exp`
#' at index 1. While you can set the assayName in these functions, it is often
#' easier to set the default assay at index 1 instead. This function aids in
#' doing so by swapping the positions of the assays when calling [assayNames].
#' @returns SummarizedExperiment object with the given assay at index 1
#' @param exp SummarizedExperiment object obtained after alignment
#' @param default Name of the assay to be set as default assay. Must be
#' present in [assayNames]
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Impute missing values
#' sumRnegative <- imputation(sumRnegative, setDefault = FALSE)
#'
#' # Set the default assay to the Imputed assay
#' sumRnegative <- setDefaultAssay(sumRnegative, "Imputed")
#'
#' # Print the assayNames
#' assayNames(sumRnegative)
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
#' @description This function sets the value of `metadata(exp)$phenotype` in
#' sumR Experiment object after validation.
#' @details The `phenotype` slot is used during post-processing, statistics,
#' and modelling to determine the groups of the [colData]. This allows for easy
#' processing by setting the phenotype column once. However, in cases of
#' multiple phenotypes present, switching the default phenotype may be a valid
#' strategy. This function aids in doing so by checking if A) the supplied
#' experiment is valid, B) the given phenotype is present as a column in
#' [colData], and C) Return the updated experiment.
#' @returns Updated SummarizedExperiment with the phenotype slot changed to
#' the supplied phenotype.
#' @param exp SummarizedExperiment object obtained after alignment
#' @param phenotype Character vector specifying the column in colData to
#' be used as phenotype
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Read the cellData
#' cellData <- read.csv(system.file("cellData.csv", package = "sumR"),
#'                      sep = ",", row.names = 1)
#'
#' # Add the cellData
#' sumRnegative <- addCellData(sumRnegative, cellData)
#'
#' # Set the phenotype
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
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

#' @title Add the Sample/Blank ratio for each compound
#' @description This function calculates the ratio between the median sample
#' area and the median blank area. The result is stored in the [rowData] slot.
#' @details The sample/blank ratio is a metric that aids in determining the
#' threshold for the blank filter. In cases with a lot of imputed values,
#' this ratio might be (very) high. In such cases, it might be suitable to
#' use a log2 or log10 assay to calculate the log fold change instead.
#' @returns A SummarizedExperiment object with results stored in the [rowData]
#' as `SampleBlankRatio`.
#' @param exp SummarizedExperiment Object
#' @param typeColumn Name of the column in `colData(exp)` that represents
#' the sample type and contains at least the values `SAMPLE` and `BLANK`.
#' @export
#' @examples
#' # Read data
#' data("sumRnegative")
#'
#' # Read the cellData
#' cellData <- read.csv(system.file("cellData.csv", package = "sumR"),
#'                      sep = ",", row.names = 1)
#'
#' # Add the cellData
#' sumRnegative <- addCellData(sumRnegative, cellData)
#'
#' # Calculate the Sample / Blank ratio
#' sumRnegative <- addSampleBlankRatio(sumRnegative)
#'
#' # Show the ratio of the first 5 compounds
#' rowData(sumRnegative[1:5, ])
addSampleBlankRatio <- function(exp, assay = 1, typeColumn = "Type"){
  samps <- as.matrix(assay(exp[, exp$Type == "SAMPLE"], assay))
  blanks <- as.matrix(assay(exp[, exp$Type == "BLANK"], assay))
  rowData(exp)$SampleBlankRatio <- round(rowMedians(samps, na.rm = TRUE) /
                                           rowMedians(blanks, na.rm = TRUE), 3)
  exp
}

