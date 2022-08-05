doBinning <- function(spectra, split = "scan", tolerance = 0.002) {
  if (class(spectra) == "data.frame") spectra <- split.data.frame(spectra, spectra[, split])

  nonEmpty <- sapply(spectra, nrow) != 0L # checking if the list is not empty
  samples <- rep.int(
    seq_along(spectra),
    sapply(spectra, nrow)
  )

  mass <- unname(unlist((lapply(
    spectra[nonEmpty],
    function(x) as.double(x$mz)
  )),
  recursive = FALSE, use.names = FALSE
  ))
  intensities <- unlist(lapply(
    spectra[nonEmpty],
    function(x) as.double(x$i)
  ),
  recursive = FALSE, use.names = FALSE
  )
  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]

  mass <- binning(
    mass = mass, intensities = intensities,
    samples = samples, tolerance = tolerance
  )

  s <- sort.int(mass, index.return = TRUE) # sort results into mass lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  lIdx <- split(seq_along(mass), samples)

  spectra[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses
    p <- NULL
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    as.data.frame(p)
  }, p = spectra[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  spectra
}

prepareSpectra <- function(spectra){
  non_nulls <- !vapply(spectra, is.null, logical(1))
  spectra <- spectra[non_nulls]
  df <- data.frame(scan = rep(which(non_nulls), sapply(spectra, nrow)), do.call(rbind, spectra))
  if (nrow(df) == 0) return(NULL)
  df
}

#' @title Bin spectra with similar mass ranges
#' @param peakList
#' @param fraction
#' @param npeaks
#' @param meanSNR
#' @param tolerance
#' @importFrom pbapply pblapply
#' @export
spectraBinning <- function(fileList, ppm = 5) {
  ppm <- ppm * 1e-6
  prepareCells(pbapply::pblapply(fileList, function(spectra){
    if (is.null(spectra)) return(NULL)
    n <- length(unique(spectra$scan))
    groups <- rep(1, n)
    df_list <- doBinning(spectra, split = "scan", tolerance = ppm)
    df <- do.call(rbind, df_list)
    rownames(df) <- 1:nrow(df)

    df$mzdiff <- df$mz - spectra$mz
    df <- df[order(df$mz), ]
    bins <- split.data.frame(df, df$mz)

    res <- DataFrame(
      mzmed = unique(df$mz),
      mzmin = sapply(bins, function(x) min(x$mz - x$mzdiff)),
      mzmax = sapply(bins, function(x) max(x$mz - x$mzdiff)),
      npeaks = sapply(bins, nrow)
    )
    res$peakidx <- lapply(bins, function(x) as.integer(rownames(x)))
    res$i <- vapply(res$peakidx, function(idx) sum(spectra$i[idx]), double(1))
    res$noise <- vapply(res$peakidx, function(idx) sum(spectra$Noise[idx]), double(1))
    res
  }))
}

spectraClustering <- function(fileList, ppm = 5){
  prepareCells(pbapply::pblapply(fileList, function(l) {
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
    res$noise <- vapply(res$peakidx, function(idx) sum(df$Noise[idx]), double(1))
    res[, -grep("rt|X1", colnames(res))]
  }))
}


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

validatePeak <- function(peakDf){
  att <- attributes(peakDf)
  all(
    class(peakDf) == "data.frame",
    all(c("scan", "mz", "i", "Noise") %in% colnames(peakDf)),
    nrow(peakDf) > 0,
    all(c("polarity", "massWindow", "combineSpectra", "Files") %in% names(att)),
    any(c("-", "+") %in% att$polarity),
    any(c(TRUE, FALSE) %in% att$combineSpectra),
    length(att$Files) == 1,
    file.exists(att$Files)
  )
}

validatePeaks <- function(peakList){
  all(
    length(names(peakList)) == length(peakList),
    all(sapply(peakList, validatePeak))
  )
}

#' @title Group spectra
#' @export
spectraAlignment <- function(peaks, method = "binning", ppm = 5, nPeaks = 0, ...){
  if (!validatePeaks(peaks)) {
    warning("Invalid peakList object")
    return(NULL)
  }
  cells <- switch(method,
         "binning" = spectraBinning(peaks, ppm = ppm,...),
         "clustering" = spectraClustering(peaks, ppm = ppm, ...),
         "density" = spectraDensity(peaks, ppm = ppm, ...)
  )
  cells[cells$npeaks >= nPeaks, ]
}

validateSpectra <- function(spectraDf){
  all(
    class(spectraDf) == "DFrame",
    attr(spectraDf, "package") == "S4Vectors",
    all(c("sample", "mz", "mzmin", "mzmax", "npeaks",
          "peakidx", "i", "noise") %in% colnames(spectraDf)),
    typeof(spectraDf$sample) == "character",
    typeof(spectraDf$mz) == "double",
    typeof(spectraDf$peakidx) == "list",
    typeof(spectraDf$i) == "double",
    nrow(spectraDf) > 0
  )
}

#' @title Group cells into SE
#' @export
cellAlignment <- function(spectraDf, method = "binning", ppm = 5, nCells = 0,
                          cellData = NULL, phenotype = NULL, ...){

  if (is.null(cellData)) {
    warning("\rNo cellData added, which is needed for postprocessing.\n",
            "Use 'addCellData()' to add cellData before continuing.")
  }

  if (!validateSpectra(spectraDf)){
    warning("Invalid aligned spectra object")
    return(NULL)
  }
  cells <- switch(method,
          "binning" = cellBinning(spectraDf, ppm = ppm, ...),
          "clustering" = cellClustering(spectraDf, ppm = ppm, ...),
          "density" = cellDensity(spectraDf, ppm = ppm, ...)
  )
  if (!is.null(phenotype)) metadata(cells)$phenotype <- phenotype

  if (!is.null(cellData)){
    idx <- intersect(colnames(cells), rownames(cellData))
    colData(cells) <- DataFrame(cellData[idx, ,drop = FALSE])
  }


  cells[rowData(cells)$npeaks >= nCells, ]
}


prepareCells <- function(cells){
  non_nulls <- !vapply(cells, is.null, logical(1))
  cells <- cells[non_nulls]

  df <- data.frame(sample = rep(names(cells), sapply(cells, nrow)), do.call(rbind, cells))
  if (nrow(df) == 0) return(NULL)
  df$rt <- -1

  colnames(df)[which(colnames(df) == "mzmed")] <- "mz"
  if (any(grepl("rt", colnames(df)))) {
    colnames(df)[which(colnames(df) == "rtmed")] <- "rt"
  } else {
    df$rt <- -1
  }
  DataFrame(df)
}

constructSE <- function(binned, spectra){

  groups <- as.integer(as.factor(rep( names(binned$npeaks), binned$npeaks)))
  spectra <- spectra[unlist(binned$peakidx), ]
  spectra$group <- groups

  df <- as.data.frame(spectra[, c("sample", "i", "group")])
  intensity <- as.matrix(DataFrame(tidyr::pivot_wider(df, names_from = "sample",
                                            values_from = "i")[,-1]))
  dimnames(intensity) <- list(1:nrow(intensity), unique(df$sample))

  df <- as.data.frame(spectra[, c("sample", "noise", "group")])
  noise <- as.matrix(DataFrame(tidyr::pivot_wider(df, names_from = "sample",
                                    values_from = "noise")[,-1]))
  dimnames(noise) <- list(1:nrow(noise), unique(df$sample))


  se <- SummarizedExperiment(
    assays = list(Area = as.matrix(intensity), Noise = as.matrix(noise)),
    rowData = DataFrame(binned),
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
cellBinning <- function(spectra, cellData = NULL, phenotype = NULL,
                        ppm = 5) {

  ppm <- ppm * 1e-6
  df_list <- sumR:::doBinning(as.data.frame(spectra), split = "sample", tolerance = ppm)

  df <- do.call(rbind, df_list)
  rownames(df) <- 1:nrow(df)

  df$mzdiff <- df$mz - spectra$mz

  df <- df[order(df$mz), ]
  bins <- split.data.frame(df, df$mz)

  res <- DataFrame(
    mz = unique(df$mz),
    mzmin = sapply(bins, function(x) min(x$mz - x$mzdiff)),
    mzmax = sapply(bins, function(x) max(x$mz - x$mzdiff)),
    npeaks = sapply(bins, nrow)
  )

  res$peakidx <- lapply(bins, function(x) as.integer(rownames(x)))

  return(constructSE(res, spectra))
}

cellClustering <- function(spectra, ppm = 5){
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

