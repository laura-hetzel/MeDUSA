#' @title Format scans in a file
#' @param files
#' @param massWindow
#' @param polarity
#' @param combineSpectra
#' @export
prepareFile <- function(file, massWindow = c(0, Inf), centroid = F,
                        polarity = "-", combineSpectra = FALSE, cl = NULL) {
  if (!is.null(cl)) {
    pbo <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(pbo), add = TRUE)
  }
  file <- file.path(file)

  if (!file.exists(file)) {
    warning(sprintf("Cannot find file %s", file))
    return(NULL)
  }
  z <- mzR::openMSfile(file)
  x <- mzR::peaks(z)
  df <- mzR::header(z)

  range <- abs(df$scanWindowUpperLimit - df$scanWindowLowerLimit)
  scans <- range >= massWindow[1] & range <= massWindow[2]

  if (!is.null(polarity)) {
    polarity_filter <- grepl("FTMS - p NSI", df$filterString)
    if (polarity == "+") polarity_filter <- !polarity_filter
    scans <- scans & polarity_filter
  }

  if (combineSpectra) {
    s <- split(which(scans), c(1, cumprod(diff(which(scans)))))
    return(lapply(s, function(z) {
      df <- do.call(rbind, x[z])
      if (centroid) {
        df <- centroid(df)
      }
      df
    }))
  } else {
    l <- x[scans]
    if (centroid) {
      l <- pbapply::pblapply(l, cl = cl, centroid)
    }
  }
  l <- prepareSpectra(l)
  if (length(l) == 0) l <- list()
  attr(l, "polarity") <- polarity
  attr(l, "massWindow") <- massWindow
  attr(l, "combineSpectra") <- combineSpectra
  attr(l, "Files") <- file
  l
}

#' @title Format scans in files
#' @param files
#' @param ... see `prepareFile`
#' @export
extractPeaks <- function(files, cores = 1, ...){

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
  result <- pbapply::pblapply(files, function(x) prepareFile(x, cl = cl, ...))
  if (!is.null(cl)) parallel::stopCluster(cl)

  # Set names of results
  result <- setNames(result, tools::file_path_sans_ext(basename(files)))

  # Return only spectra with found centroid apexes
  result[!vapply(result, function(x) is.null(x) | length(x) == 0, logical(1))]
}

#' @title Format scans in a directory
#' @param directory
#' @param ... see `prepareFile`
#' @export
extractPeaksBatch <- function(directory, ...){

  # Check if directory exists
  directory <- dirname(directory)
  if (!dir.exists(directory)) {
    warning(sprintf("Cannot find directory %s", directory))
    return(NULL)
  }

  # Check if mzML files are found in directory
  files <- list.files(directory, pattern = ".mzML", full.names = T, recursive = F)
  if (length(files) == 0) {
    warning(sprintf("Cannot find mzML files in directory %s", directory))
    return(NULL)
  }

  # Prepare mzMLs
  extractPeaks(files, ...)
}



calculateNoise <- function(spectrum, centroids, noiseWindow){
  mzInd <- spectrum[, 1]
  noise <- cwt_new(spectrum[, 2], scales = 1)
  minNoiseLevel <- max(noise) * 0.001
  noise <- abs(noise)
  nMz <- length(mzInd)
  winSize.noise <- as.integer(nMz * noiseWindow)
  sapply(centroids, function(k){
    ind.k <- mzInd[k]
    start.k <- max(1, ind.k - winSize.noise)
    end.k <- min(nMz, ind.k + winSize.noise)

    noiseLevel.k <- quantile(noise[start.k:end.k], probs = 0.95, na.rm = T)
    if (noiseLevel.k < minNoiseLevel) {
      noiseLevel.k <- minNoiseLevel
    }
    noiseLevel.k
  })

}

centroid <- function(spectrum, halfWindowSize = 2L, noiseWindow = 0.0001) {
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
savgol <- function(y, halfWindowSize = 10L, polynomialOrder = 3L) {
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
