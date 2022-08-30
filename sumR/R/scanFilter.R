#' @title Perform centroiding of scans in peaks
#' @description This function converts profile-mode data into centroided
#' data. In case of the use of different mass-windows, the `scans` parameter
#' can be set to subset the scans to be used. The parameter `combineSpectra`
#' can be set to TRUE if you want to combine the multiple mass-windows into
#' a single full-range scan.
#' @details Centroiding is a major part of a typical metabolomics pipeline.
#' It reduces the data from a profile-mode dataset drastically, by picking
#' the local maxima of profiles after smoothing. While in high-resolution
#' datasets peak picking is followed afterwards, for single-cell DI-MS, peak
#' picking after centroiding provides sub-optimal results. In many cases, the
#' peak picking would find peaks originating from the solution, instead of the
#' compounds in the cells. Therefore, it is recommended to only perform
#' centroiding and remove the solution-originating peaks afterwards using
#' filtering.
#' @param peaks List of dataframes, usually obtained from `peaks()` of the
#' `mzR` package.
#' @param scans Vector of indices, which scans should be used for centroiding?
#' Defaults to all scans in `peaks`
#' @param combineSpectra Boolean value, should scans of different masswindows
#' be combined into a full-range scan before centroiding? Defaults to `FALSE`
#' @returns  data.frame of centroided peaks with the columns `scan`, `mz`, `i`
#' and `Noise`
#' @importFrom pbapply pblapply
#' @export
doCentroid <- function(peaks, scans = seq_len(length(x)),
                       combineSpectra = FALSE, cl = NULL){
  if (combineSpectra) {
    s <- split(which(scans), c(1, cumprod(diff(which(scans)))))
    l <- lapply(s, function(z) {
      centroid(do.call(rbind, peaks[z]))
    })
  } else {
    l <- pbapply::pblapply(peaks[scans], cl = cl, centroid)
  }
  formatSpectra(l)
}

#' @title Filter non-empty spectra and format into a dataframe
#' @description This function removes any spectra without centroided peaks
#' and combines the spectra from a list into a dataframe with an added
#' column called `scans`
#' @param spectra List of centroided spectra, usually obtained from
#' `doCentroid`
#' @returns data.frame with centroided peaks.
#' @noRd
formatSpectra <- function(spectra){
  non_nulls <- !vapply(spectra, is.null, logical(1))
  spectra <- spectra[non_nulls]
  scans <- rep(which(non_nulls), vapply(spectra, nrow, integer(1)))
  df <- data.frame(scan = scans, do.call(rbind, spectra))
  if (nrow(df) == 0) return(NULL)
  df
}

#' @title Centroid a provided mzML file
#' @description Centroiding is a vital part of a typical metabolomics workflow.
#' This function centroids a given mzML file and returns a dataframe of
#' centroided peaks. See `doCentroid` for more information about centroiding.
#' @param file A mzML file in profile-mode
#' @param massWindow A vector of size 2 with min-max size of mass-windows to be
#' used. E.g. A mass-window from 200 to 300 has a mass-window of 100. Defaults
#' to c(0, Inf) indicating that all mass-windows should be taken into account
#' @param polarity A single character of either "-" (negative) or "+"
#' (positive) indicating the polarity to be used. Defaults to "-"
#' @param combineSpectra Boolean value, should scans of different masswindows
#' be combined into a full-range scan before centroiding? Defaults to `FALSE`
#' @returns  data.frame of centroided peaks with the columns `scan`, `mz`, `i`
#' and `Noise`. It also stores the used parameters as attributes of the
#' data.frame
#' @importFrom pbapply pboptions
#' @importFrom mzR openMSfile header peaks
#' @inherit doCentroid details
#' @export
centroidFile <- function(file, massWindow = c(0, Inf), polarity = "-",
                         combineSpectra = FALSE, cl = NULL) {
  if (!is.null(cl)) {
    pbo <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(pbo), add = TRUE)
  }
  file <- file.path(file)

  if (!file.exists(file)) {
    warning(sprintf("Cannot find file %s", file))
    return(NULL)
  }
  z <- openMSfile(file)
  df <- header(z)

  range <- abs(df$scanWindowUpperLimit - df$scanWindowLowerLimit)
  scans <- range >= massWindow[1] & range <= massWindow[2]

  if (!is.null(polarity)) {
    polarity_filter <- grepl("FTMS - p NSI", df$filterString)
    if (polarity == "+") polarity_filter <- !polarity_filter
    scans <- scans & polarity_filter
  }

  l <- doCentroid(peaks(z), scans, combineSpectra = combineSpectra, cl = cl)
  if (length(l) == 0) l <- list()
  attr(l, "polarity") <- polarity
  attr(l, "massWindow") <- massWindow
  attr(l, "combineSpectra") <- combineSpectra
  attr(l, "Files") <- file
  l
}

#' @title Format scans in files
#' @param files A vector of paths to mzML files.
#' @param cores Integer value of the number of cores to use. Any value > 1
#' indicates a multi-core approach. Defaults to 1.
#' @inheritDotParams centroidFile -file
#' @inherit centroidFile details
#' @importFrom parallel makeCluster stopCluster
#' @importFrom pbapply pblapply
#' @importFrom tools file_path_sans_ext
#' @returns A list of data.frames, each obtained from `centroidFile`. The names
#' of the list are set to the basename of each file provided in `files`
#' @export
#' @examples
#' # Define files
#' folder <- system.file("cells", package = "sumR")
#' files <- list.files(folder, full.names = TRUE)[seq_len(4)]
#'
#' # Extract Peaks
#' extractPeaks(files, polarity = "-")
#'
#' # Extract Peaks in parallel
#' extractPeaks(files, polarity = "-", cores = 4)
#'
#' # Extract Peaks of positive polarity
#' extractPeaks(files, polarity = "+")
extractPeaks <- function(files, cores = 1, ...){
  files <- file.path(files)
  # Check if mzML files are valid files
  files <- files[file.exists(files)]
  files <- files[grep(".mzML", files)]

  if (length(files) == 0) {
    warning("Cannot find mzML files in given files")
    return(NULL)
  }

  # Make cluster if cores > 1
  cl <- NULL
  if (cores > 1) {
    cl <- makeCluster(cores)
  }

  # Prepare each file
  result <- pblapply(files, function(x) centroidFile(x, cl = cl, ...))
  if (!is.null(cl)) stopCluster(cl)

  # Set names of results
  result <- setNames(result, file_path_sans_ext(basename(files)))

  # Return only spectra with found centroid apexes
  result[!vapply(result, function(x) is.null(x) | length(x) == 0, logical(1))]
}

#' @title Format scans in a directory
#' @param directory Path to directory containing .mzML files in profile-mode.
#' @inheritDotParams extractPeaks -files
#' @inheritDotParams centroidFile -file
#' @inherit centroidFile details
#' @inherit extractPeaks return
#' @export
extractPeaksBatch <- function(directory, ...){

  # Check if directory exists
  directory <- dirname(directory)
  if (!dir.exists(directory)) {
    warning(sprintf("Cannot find directory %s", directory))
    return(NULL)
  }

  # Check if mzML files are found in directory
  files <- list.files(directory, pattern = ".mzML",
                      full.names = T, recursive = F)
  if (length(files) == 0) {
    warning(sprintf("Cannot find mzML files in directory %s", directory))
    return(NULL)
  }

  # Prepare mzMLs
  extractPeaks(files, ...)
}


#' @title Calculate the noise in a given spectrum
#' @param spectrum
#' @param centroids
#' @param noiseWindow
#' @importFrom stats quantile
#' @noRd
calculateNoise <- function(spectrum, centroids, noiseWindow){
  mzInd <- spectrum[, 1]
  noise <- cwt_new(spectrum[, 2], scales = 1)
  minNoiseLevel <- max(noise) * 0.001
  noise <- abs(noise)
  nMz <- length(mzInd)
  winSize.noise <- as.integer(nMz * noiseWindow)
  vapply(centroids, function(k){
    ind.k <- mzInd[k]
    start.k <- max(1, ind.k - winSize.noise)
    end.k <- min(nMz, ind.k + winSize.noise)

    noiseLevel.k <- quantile(noise[start.k:end.k], probs = 0.95, na.rm = T)
    if (noiseLevel.k < minNoiseLevel) {
      noiseLevel.k <- minNoiseLevel
    }
    noiseLevel.k
  }, double(1))
}

#' @title Centroid a given spectrum
#' @param spectrum
#' @param halfWindowSize
#' @param noiseWindow
#' @noRd
centroid <- function(spectrum, halfWindowSize = 2L, noiseWindow = 0.0001) {
  if (nrow(spectrum) == 0) return(NULL)
  intensities <- savgol(spectrum[, 2], halfWindowSize)
  intensities[intensities < 0] <- 0

  ## find local maxima
  centroids <- which(diff(sign(diff(intensities))) == -2) + 1

  peakNoise <- calculateNoise(spectrum, centroids, noiseWindow)
  spectrum <- as.data.frame(spectrum[centroids, ])
  colnames(spectrum) <- c("mz", "i")
  spectrum$Noise <- peakNoise
  spectrum
}

#' @title Apply Savitz-golay filter
#' @param y Vector of intensities
#' @param halfWindowSize
#' @param polynomialOrder
#' @noRd
savgol <- function(y, halfWindowSize = 2L, polynomialOrder = 3L) {

  sav.filter <- function(x, hws, coef) {
    n <- length(x)
    w <- 2L * hws + 1L
    y <- stats::filter(x = x, filter = coef[hws + 1L, ], sides = 2L)
    attributes(y) <- NULL
    y[seq_len(hws)] <- head(coef, hws) %*% head(x, w)
    y[(n - hws + 1L):n] <- tail(coef, hws) %*% tail(x, w)
    y
  }

  coef <- function(m, k = 3L) {
    k <- 0L:k
    nm <- 2L * m + 1L
    nk <- length(k)
    K <- matrix(k, nrow = nm, ncol = nk, byrow = TRUE)
    F <- matrix(double(), nrow = nm, ncol = nm)
    for (i in seq_len(m + 1L)) {
      M <- matrix(seq_len(nm) - i, nrow = nm, ncol = nk, byrow = FALSE)
      X <- M^K
      T <- solve(t(X) %*% X) %*% t(X)
      F[i, ] <- T[1L, ]
    }
    F[(m + 2L):nm, ] <- rev(F[seq_len(m), ])
    F
  }

  sav.filter(y, hws = halfWindowSize, coef = coef(m = halfWindowSize,
                                                  k = polynomialOrder))
}
