#' @title Bin masses of spectra
#' @description This is an adapted form of Maldiquants binning function.
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
#' @returns A list of spectra with the masses replaced by the mean of the
#' associated bin
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

#' @title Bin spectra with similar mass ranges
#' @description
#' @returns DataFrame object from S4Vectors containing aligned spectra
#' @param spectraList List of spectra to be binned
#' @param ppm How much can masses deviate from the mean of the bin?
#' Defaults to 5 ppm
#' @inherit doBinning details
#' @importFrom pbapply pblapply
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
    res
  })

  res <- prepareCells(res)
  attr(res, "Datetime") <- lapply(spectraList, function(x) attr(x, "Datetime"))
  res
}

#' @title Group spectra using clustering
#' @param spectraList List of spectra to be binned
#' @param ppm How much can masses deviate from the mean of the cluster?
#' Defaults to 5 ppm
#' @inherit spectraBinning details
#' @importFrom pbapply pblapply
#' @importFrom S4Vectors DataFrame
#' @export
spectraClustering <- function(spectraList, ppm = 5){
  if (requireNamespace("xcms", quietly = TRUE)) {
    prepareCells(pbapply::pblapply(spectraList, function(l) {
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

#' @title Group spectra
#' @export
spectraAlignment <- function(peaks, method = "binning", ppm = 5, nPeaks = 0){
  if (!validatePeaks(peaks)) {
    warning("Invalid peakList object")
    return(NULL)
  }
  cells <- switch(method,
         "binning" = spectraBinning(peaks, ppm = ppm),
         "clustering" = spectraClustering(peaks, ppm = ppm)
  )
  spectra <- cells[cells$npeaks >= nPeaks, ]
  spectra <- spectra[order(spectra$mz), ]
  rownames(spectra) <- 1:nrow(spectra)
  spectra
}



#' @title Group cells into SE
#' @importFrom lubridate as_datetime
#' @export
cellAlignment <- function(spectraDf, method = "binning", ppm = 5, nCells = 0,
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

  datetime <- data.frame(Datetime = do.call(rbind, attr(spectra, "Datetime")))
  datetime <- datetime[order(rownames(datetime)), , drop = FALSE]
  cells <- cells[, order(colnames(cells))]

  colData(cells) <- cbind(colData(cells), datetime)
  colData(cells)$Datetime <- lubridate::as_datetime(colData(cells)$Datetime)
  cells <- cells[, order(colData(cells)$Datetime)]
  cells <- cells[rowData(cells)$npeaks >= nCells, ]
  rownames(cells) <- 1:nrow(cells)
  cells
}


prepareCells <- function(cells){
  non_nulls <- !vapply(cells, is.null, logical(1))
  cells <- cells[non_nulls]

  df <- data.frame(sample = rep(names(cells), sapply(cells, nrow)), do.call(rbind, cells))
  if (nrow(df) == 0) return(NULL)
  DataFrame(df)
}

constructSE <- function(res, spectra){

  spectra$group <- as.integer(as.factor(rep( names(res$npeaks), res$npeaks)))
  df <- as.data.frame(spectra[, c("sample", "i", "group")])
  intensity <- as.matrix(DataFrame(tidyr::pivot_wider(df, names_from = "sample",
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

#' @title Bin cells
#' @param spectraList List of spectra
#' @param tolerance
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
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
                                  sampleGroups = groups, ppm = ppm,
                                  minFraction = 0)
    )

    res <- DataFrame(peaks$featureDefinitions)
    res$peakidx <- peaks$peakIndex

    res <- res[, grep("mz|npeaks|peakidx", colnames(res))]
    return(constructSE(DataFrame(res), spectra))
  }
}

