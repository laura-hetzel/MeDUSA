#' @title Format scans in a file
#' @param files
#' @param massWindow
#' @param polarity
#' @param combineSpectra
#' @export
prepareFile <- function(file, massWindow = c(0, Inf), centroid = F,
                        polarity = "-", combineSpectra = FALSE) {

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
  }
  l <- x[scans]
  if (centroid) l <- lapply(l, centroid)
  l
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
prepareFiles <- function(files, a = 2, ...){
  result <- pbapply::pblapply(files, function(x) prepareFile(x, ...))
  setNames(result, tools::file_path_sans_ext(basename(files)))
}

#' @title Format scans in a directory
#' @param directory
#' @param ... see `prepareFile`
#' @export
prepareBatch <- function(directory, ...){
  files <- list.files(directory, pattern = ".mzML", full.names = T, recursive = F)
  prepareFiles(files, ...)
}





centroid <- function(spectrum, halfWindowSize = 2L) {
  intensities <- savgol(spectrum[, 2], halfWindowSize)
  intensities[intensities < 0] <- 0

  ## find local maxima
  spectrum[which(diff(diff(intensities) >= 0) < 0) + 1, ]
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
