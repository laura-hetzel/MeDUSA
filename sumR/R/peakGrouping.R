#' @title Group signals across spectra into peaks
#' @description This function will perform either binning (default) or
#' clustering to group signals across spectra into peaks from a single file
#' @details With DI-MS, peak in retention time space are too unreliable to
#' properly peak pick. Here, we use binning to bin mass signals that deviate
#' less than the given ppm.
#' @returns DataFrame where each row is a peak
#' @param peaks Centroided spectra obtained from [extractPeaks].
#' @param method Character, one of either `"binning"` (default) or
#' `"clustering"`.
#' @param ppm Tolerance between masses to be considered part of the same bin
#' or cluster in parts per million (ppm). Defaults to 5.
#' @export
#' @examples
#' # Define files
#' folder <- system.file("cells", package = "sumR")
#' files <- list.files(folder, full.names = TRUE)
#'
#' # Extract peaks from files
#' peaks <- extractPeaks(files)
#'
#' # Align peaks across spectra
#' spectra <- spectraAlignment(peaks)
#'
#' # Show first few peaks
#' print(head(spectra))
spectraAlignment <- function(peaks, method = "binning", ppm = 5){
  if (!validatePeaks(peaks)) {
    warning("Invalid peakList object")
    return(NULL)
  }
  spectra <- switch(method,
                    "binning" = spectraBinning(peaks, ppm = ppm),
                    "clustering" = spectraClustering(peaks, ppm = ppm)
  )
  spectra <- spectra[order(spectra$mz), ]
  rownames(spectra) <- 1:nrow(spectra)
  spectra
}

#' @title Bin masses of spectra
#' @description This is an adapted form of Maldiquants binning function. It will
#' check the given masses and if subsequent masses fall within the ppm
#' tolerance, it will assign them to the same bin.
#' @details Binning is a form of alignment that groups together masses
#' within a given tolerance. Due to measurements errors in mass spectrometry,
#' binning is necessary to determine which mass intensities belong to
#' the same peak.
#' @param spectra a data.frame or list of data.frames containing centroided
#' spectra.
#' @param split In case of a data.frame in `spectra`, which column should
#' be used to indicate different scans on which the alignment should occur?
#' Defaults to "scan"
#' @param ppm How much can masses deviate from the mean of the bin?
#' Defaults to 5 ppm
#' @returns A data.frame with the binned masses, their intensities and
#' original peak identifiers for traceback.
doBinning <- function(spectra, split = "scan", ppm = 5) {
  tolerance <- ppm * 1e-6
  if (is(spectra, "data.frame")) {
    spectra <- split.data.frame(spectra, spectra[, split])
  }

  # checking if the list is not empty
  nonEmpty <- as.vector(vapply(spectra, nrow, integer(1)) != 0L)
  samples <- rep.int(seq_along(spectra), vapply(spectra, nrow, double(1)))
  specs <- spectra[nonEmpty]

  mass <- unname(unlist(lapply(specs, function(x) as.double(x$mz))))
  intensities <- unlist(lapply(specs, function(x) as.double(x$i)))
  rts <- unlist(lapply(specs, function(x) as.double(x$rt)))

  # sort vector based on masses lowest to highest
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  rts <- rts[s$ix]

  ids <- binning(
    mass = mass, intensities = intensities,
    samples = samples, tolerance = tolerance
  )

  data.frame(
    mass = mass,
    rt = rts,
    intensities = intensities,
    samples = samples,
    peakIdx = as.integer(as.factor(ids)),
    order = s$ix
  )
}

#' @title Bin spectra with similar masses
#' @inherit doBinning description details
#' @returns DataFrame object from S4Vectors containing aligned spectra
#' @param spectraList List of spectra to be binned
#' @param ppm How much can masses deviate from the mean of the bin?
#' Defaults to 5 ppm
#' @importFrom pbapply pblapply
#' @importFrom stats median
#' @export
spectraBinning <- function(spectraList, ppm = 5) {
  res <- pbapply::pblapply(spectraList, function(spectra){
    if (is.null(spectra)) return(NULL)
    n <- length(unique(spectra$scan))

    df <- doBinning(spectra, split = "scan", ppm = ppm)
    df$origId <- rownames(spectra)[df$order]
    bins <- split.data.frame(df, df$peakIdx)

    res <- DataFrame(
      mz = sapply(bins, function(x) median(x$mass)),
      mzmin = sapply(bins, function(x) min(x$mass)),
      mzmax = sapply(bins, function(x) max(x$mass)),
      rt = sapply(bins, function(x) median(x$rt)),
      rtmin = sapply(bins, function(x) min(x$rt)),
      rtmax = sapply(bins, function(x) max(x$rt)),
      npeaks = sapply(bins, nrow),
      i = sapply(bins, function(x) sum(x$intensities))
    )
    res$fpeaks <- round(res$npeaks / n * length(attr(spectra, "scanranges")), 3)
    res$ppm <- round((res$mzmax - res$mzmin) / res$mz * 1e6, 3)
    res$peakIdx <- lapply(bins, function(x) x$origId)
    res[res$i > 0, ]
  })

  res <- prepareCells(res)
  attr(res, "Datetime") <- lapply(spectraList, function(x) attr(x, "Datetime"))
  res
}

#' @title Group spectra using clustering
#' @inheritParams spectraAlignment
#' @inherit spectraBinning details
#' @importFrom pbapply pblapply
#' @importFrom S4Vectors DataFrame
spectraClustering <- function(peaks, ppm = 5){
  if (requireNamespace("xcms", quietly = TRUE)) {
    prepareCells(pbapply::pblapply(peaks, function(l) {
      df <- prepareSpectra(l)
      if (is.null(df)) return(NULL)

      n <- length(unique(df$sample))
      groups <- rep(1, n)

      peaks <- quiet(xcms::do_groupPeaks_mzClust(df, sampleGroups = groups,
                                                 ppm = ppm, minFraction = 0)
      )
      res <- DataFrame(peaks$featureDefinitions)
      res$peakidx <- peaks$peakIndex

      res$i <- vapply(res$peakidx, function(idx) sum(df$i[idx]), double(1))
      res[, -grep("X1", colnames(res))]
    }))
  }
}





#' @title Align spectra across cells using binning or clustering
#' @description This function will perform either binning (default) or
#' clustering to group peaks from spectra across cells.
#' @details Alignment across cells is a vital part of metabolomics data analysis
#' in order to make comparisons across cells. Here, we support two methods of
#' alignment: binning and clustering. Binning is an adjustment from Maldiquant
#' while clustering is using xcms instead. Both will result in a
#' SummarizedExperiment object that can be used for postprocessing.
#' @returns SummarizedExperiment object that contains both data and metadata.
#' @param spectraDf DataFrame of peaks in spectra, where 1 row equals 1 peak
#' from a single cell
#' @param method Character, one of either `"binning"` (default) or
#' `"clustering"`.
#' @param ppm Tolerance between masses to be considered part of the same bin
#' or cluster in parts per million (ppm). Defaults to 5.
#' @param cellData data.frame of metadata of each cell. While optional in
#' this step, it is crucial for post-processing. Defaults to `NULL`
#' @param phenotype Name of the column in `cellData` that describes the
#' phenotype of each cell. Mandatory for statistics and modelling, but can
#' be left empty here. Defaults to `NULL`
#' @importFrom lubridate as_datetime
#' @export
#' @examples
#' # Define files
#' folder <- system.file("cells", package = "sumR")
#' files <- list.files(folder, full.names = TRUE)
#'
#' # Extract peaks from files
#' peaks <- extractPeaks(files)
#'
#' # Align peaks across spectra
#' spectra <- spectraAlignment(peaks)
#'
#' # Align spectra across cells
#' cells <- cellAlignment(spectra, method = "binning")
#'
#' ## Same but with clustering
#' # cells <- cellAlignment(spectra, method = "clustering")
#'
#' # Show output
#' print(cells)
cellAlignment <- function(spectraDf, method = "binning", ppm = 5,
                          cellData = NULL, phenotype = NULL){

  if (is.null(cellData)) {
    warning("\rNo cellData added, which is needed for postprocessing.\n",
            "Use 'addCellData()' to add cellData before continuing.")
  }

  if (!validateSpectra(spectraDf)) {
    warning("Invalid aligned spectra object")
    return(NULL)
  }
  cells <- switch(method,
          "binning" = cellBinning(spectraDf, ppm),
          "clustering" = cellClustering(spectraDf, ppm),
  )
  if (!is.null(phenotype)) metadata(cells)$phenotype <- phenotype

  if (!is.null(cellData)) cells <- addCellData(cells, cellData)

  datetime <- data.frame(Datetime = do.call(rbind, attr(spectraDf, "Datetime")))
  datetime <- datetime[order(rownames(datetime)), , drop = FALSE]
  cells <- cells[, order(colnames(cells))]

  colData(cells) <- cbind(colData(cells), datetime)
  colData(cells)$Datetime <- lubridate::as_datetime(colData(cells)$Datetime)
  cells <- cells[, order(colData(cells)$Datetime)]
  rownames(cells) <- 1:nrow(cells)
  cells
}

#' @title Combine spectra into a single DataFrame object
#' @description This function combines a list of spectra into a S4Vectors
#' DataFrame object.
#' @details Storing the spectra as a list of S4Vectors allows to store the
#' peakIdx as a list with other spectra information.
#' @returns A DataFrame object of spectra of multiple cells
#' @param cells List of spectra to be combined into a DataFrame
#' @importFrom S4Vectors DataFrame
prepareCells <- function(cells){
  non_nulls <- !vapply(cells, is.null, logical(1))
  cells <- cells[non_nulls]

  df <- data.frame(sample = rep(names(cells), sapply(cells, nrow)),
                   do.call(rbind, cells))
  if (nrow(df) == 0) return(NULL)
  DataFrame(df)
}

#' @title Construct a SummarizedExperiment from aligned spectra
#' @description This function constructs a SummarizedExperiment from a DataFrame
#' of aligned spectra.
#' @details The SummarizedExperiment is a community-standard container that
#' allows for easy requesting and manipulating data and metadata.
#' @returns SummarizedExperiment object with assay `"Area"` and rowData with
#' mz, rt and other metadata.
#' @param res DataFrame result from either [cellBinning] or [cellClustering].
#' @param spectra Large DataFrame where each row is a peak from a single sample.
#' @importFrom S4Vectors DataFrame
#' @importFrom tidyr pivot_wider
#' @importFrom SummarizedExperiment SummarizedExperiment
constructSE <- function(res, spectra){

  spectra$group <- as.integer(as.factor(rep( names(res$npeaks), res$npeaks)))
  df <- as.data.frame(spectra[, c("sample", "i", "group")])
  intensity <- as.matrix(DataFrame(pivot_wider(df, names_from = "sample",
                                                      values_from = "i")[,-1]))

  dimnames(intensity) <- list(1:nrow(intensity), unique(df$sample))
  res <- res[, -which(colnames(res) == "i"), drop = FALSE]

  se <- SummarizedExperiment(
    assays = list(Area = as.matrix(intensity)),
    rowData = DataFrame(res),
    colData = DataFrame(row.names = unique(df$sample))
  )
  rownames(se) <- 1:nrow(se)
  se
}

#' @title Bin peaks from spectra across cells.
#' @description This function will use the same binning method as
#' [spectraBinning] but across cells instead of spectra.
#' @inherit doBinning details
#' @returns A SummarizedExperiment with an Area assay and mz, rt and ppm
#' metrics in rowData.
#' @param spectra DataFrame of spectra obtained from [spectraAlignment].
#' @param ppm Tolerance between masses to be considered part of the same bin in
#' parts per million (ppm). Defaults to 5.
#' @importFrom S4Vectors DataFrame
cellBinning <- function(spectra, ppm = 5) {
  df <- doBinning(as.data.frame(spectra), split = "sample", ppm = ppm)
  rownames(df) <- 1:nrow(df)
  bins <- split.data.frame(df, df$peakIdx)

  res <- DataFrame(
    mz = sapply(bins, function(x) median(x$mass)),
    mzmin = sapply(bins, function(x) min(x$mass)),
    mzmax = sapply(bins, function(x) max(x$mass)),
    rt = sapply(bins, function(x) median(x$rt)),
    rtmin = sapply(bins, function(x) min(x$rt)),
    rtmax = sapply(bins, function(x) max(x$rt)),
    npeaks = sapply(bins, nrow),
    i = sapply(bins, function(x) sum(x$intensities))
  )

  res$peakIdx <- lapply(bins, rownames)
  res$sampleIdx <- lapply(bins, `[[`, "samples")
  res$ppm <- round((res$mzmax - res$mzmin) / res$mz * 1e6, 3)
  res$fsample <- round(res$npeaks / length(unique(spectra$sample)), 3)
  return(constructSE(res, spectra))
}

cellClustering <- function(spectra, ppm = 5){
  if (requireNamespace("xcms", quietly = TRUE)) {

    n <- length(unique(spectra$sample))
    groups <- rep(1, n)
    peaks <-  quiet(
      xcms::do_groupPeaks_mzClust(as.data.frame(spectra),
                                  sampleGroups = groups,
                                  ppm = ppm,
                                  minFraction = 0)
    )

    res <- DataFrame(peaks$featureDefinitions)
    res$peakidx <- peaks$peakIndex

    res <- res[, grep("mz|npeaks|peakidx", colnames(res))]
    return(constructSE(DataFrame(res), spectra))
  }
}

