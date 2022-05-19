#' @title Pick peaks in a spectrum
#' @importFrom MassSpecWavelet getLocalMaximumCWT getRidge
pickSpectra <- function(x, SNR = 0) {
  ms <- as.data.frame(x)
  colnames(ms) <- c("mz", "i")
  wCoefs <- cbind("0" = ms$i, cwt_new(ms$i))
  localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs)
  ridgeList <- MassSpecWavelet::getRidge(localMax)
  majorPeakInfo <- identifyPeaks(ms$i, ridgeList, wCoefs,
                                 SNR.Th = SNR, nearbyPeak = TRUE
  )

  df <- cbind(ms[as.vector(majorPeakInfo$allPeakIndex), ],
              Noise = majorPeakInfo$peakNoise[names(majorPeakInfo$allPeakIndex)],
              Value = majorPeakInfo$peakValue[names(majorPeakInfo$allPeakIndex)],
              SNR = majorPeakInfo$peakSNR[names(majorPeakInfo$allPeakIndex)]
  )
  colnames(df) <- make.names(colnames(df), unique = T)
  df[df$SNR >= SNR, ]
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

  sav.filter(y, hws = halfWindowSize, coef = coef(m = halfWindowSize, k = polynomialOrder))
}

#' @title Perform peak picking
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom mzR openMSfile peaks header
#' @importFrom dplyr distinct
#' @importFrom tools file_path_sans_ext
#' @export
peakPicking <- function(files, doCentroid = F, massDefect = 0.8, polarity = "-",
                        cores = detectCores(logical = F), SNR = 0) {
  cl <- makeCluster(cores)
  clusterExport(cl, varlist = names(sys.frame()))
  result <- pbapply::pblapply(files, cl = cl, function(f) {
    z <- mzR::openMSfile(f)
    x <- mzR::peaks(z)
    if (doCentroid) {
      df <- mzR::header(z)
      scans <- df$scanWindowUpperLimit - df$scanWindowLowerLimit <= 200
      polarity_filter <- grepl("FTMS - p NSI", df$filterString)
      if (polarity == "+") polarity_filter <- !polarity_filter
      scans <- scans & polarity_filter
      x <- lapply(setNames(x[scans], df[scans, ]$retentionTime), centroid)
    }

    l <- lapply(x, function(spectrum) {
      tryCatch(suppressWarnings(pickSpectra(spectrum, SNR)), error = function(err) NULL)
    })

    non_nulls <- !vapply(l, is.null, logical(1))
    l <- l[non_nulls]
    df <- data.frame(scan = rep(which(non_nulls), sapply(l, nrow)), do.call(rbind, l))
    df <- dplyr::distinct(df[df$Value > 0, ])
    if (nrow(df) == 0) return(NULL)
    rownames(df) <- 1:nrow(df)
    df <- MassDefectFilter(df, mz_MD_plot = F)
    df[df$MD < massDefect, ]
  })
  samps <- tools::file_path_sans_ext(basename(files))
  stopCluster(cl)
  names(result) <- samps
  result[!vapply(result, is.null, logical(1))]
}

identifyPeaks <- function(ms, ridgeList, wCoefs, scales = as.numeric(colnames(wCoefs)),
                          SNR.Th = 3, peakScaleRange = 5, ridgeLength = 32, nearbyPeak = FALSE,
                          nearbyWinSize = 100, winSize.noise = 500, SNR.method = "quantile",
                          minNoiseLevel = 0.001) {
  if (ridgeLength > max(scales)) ridgeLength <- max(scales)
  peakScaleRange <- scales[scales >= peakScaleRange]

  minNoiseLevel <- max(wCoefs) * minNoiseLevel

  ridgeLen <- sapply(ridgeList, length)
  ridgeName <- names(ridgeList)
  ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName, "_"))), nrow = 2)
  ridgeLevel <- ridgeInfo[1, ]
  notnull <- sapply(ridgeList, function(x) {
    !is.null(x[1])
  })
  mzInd <- sapply(ridgeList[notnull], function(x) {
    x[1]
  })
  ord <- order(mzInd)
  ridgeName <- ridgeName[ord]
  ridgeLen <- ridgeLen[ord]
  ridgeLevel <- ridgeLevel[ord]
  ridgeList <- ridgeList[ord]
  mzInd <- mzInd[ord]
  peakScale <- NULL
  peakCenterInd <- NULL
  peakValue <- NULL
  for (i in 1:length(ridgeList)) {
    ridge.i <- ridgeList[[i]]
    level.i <- ridgeLevel[i]
    levels.i <- level.i:(level.i + ridgeLen[i] - 1)
    scales.i <- scales[levels.i]
    selInd.i <- which(scales.i %in% peakScaleRange)
    if (length(selInd.i) == 0) {
      peakScale <- c(peakScale, scales.i[1])
      peakCenterInd <- c(peakCenterInd, ridge.i[1])
      peakValue <- c(peakValue, 0)
      next
    }
    levels.i <- levels.i[selInd.i]
    scales.i <- scales.i[selInd.i]
    ridge.i <- ridge.i[selInd.i]
    if (scales.i[1] == 0) {
      ind.i <- cbind(ridge.i[-1], levels.i[-1])
    } else {
      ind.i <- cbind(ridge.i, levels.i)
    }
    ridgeValue.i <- wCoefs[ind.i]
    maxInd.i <- which.max(ridgeValue.i)
    peakScale <- c(peakScale, scales.i[maxInd.i])
    peakCenterInd <- c(peakCenterInd, ridge.i[maxInd.i])
    peakValue <- c(peakValue, ridgeValue.i[maxInd.i])
  }
  noise <- abs(wCoefs[, "1"])
  peakSNR <- NULL
  peakNoise <- NULL
  nMz <- nrow(wCoefs)
  for (k in 1:length(ridgeList)) {
    ind.k <- mzInd[k]
    start.k <- ifelse(ind.k - winSize.noise < 1, 1, ind.k -
                        winSize.noise)
    end.k <- ifelse(ind.k + winSize.noise > nMz, nMz, ind.k +
                      winSize.noise)

    noiseLevel.k <- quantile(noise[start.k:end.k], probs = 0.95)
    if (noiseLevel.k < minNoiseLevel) {
      noiseLevel.k <- minNoiseLevel
    }
    peakSNR <- c(peakSNR, peakValue[k] / noiseLevel.k)
    peakNoise <- c(peakNoise, noiseLevel.k)
  }
  selInd1 <- (scales[ridgeLevel + ridgeLen - 1] >= ridgeLength)
  if (nearbyPeak) {
    selInd1 <- which(selInd1)
    index <- 1:length(mzInd)
    nearbyWinSize <- 150
    tempInd <- NULL
    for (ind.i in selInd1) {
      tempInd <- c(tempInd, index[mzInd >= mzInd[ind.i] -
                                    nearbyWinSize & mzInd <= mzInd[ind.i] + nearbyWinSize])
    }
    selInd1 <- (index %in% tempInd)
  }
  selInd2 <- (peakSNR > SNR.Th)
  selInd3 <- !(mzInd %in% c(1:(nearbyWinSize / 2), (nrow(wCoefs) -
                                                      (nearbyWinSize / 2) + 1):nrow(wCoefs)))
  selInd <- (selInd1 & selInd2 & selInd3)
  names(peakSNR) <- names(peakScale) <- names(peakNoise) <- names(peakCenterInd) <- names(peakValue) <- names(mzInd) <- ridgeName
  return(list(
    peakIndex = mzInd[selInd], peakValue = peakValue,
    peakCenterIndex = peakCenterInd, peakSNR = peakSNR,
    peakNoise = peakNoise,
    peakScale = peakScale,
    potentialPeakIndex = mzInd[selInd1 & selInd3], allPeakIndex = mzInd
  ))
}


cwt_new <- function(ms, scales = NULL, max_scale = 32) {
  psi_xval <- seq(-8, 8, length = 1024)
  psi <- (2 / sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2 / 2)

  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax <- psi_xval[length(psi_xval)]

  if (is.null(scales)) {
    scales <- seq(1, ifelse(floor(length(ms) / 16) < max_scale,
                            floor(length(ms) / 16), max_scale
    ), 2)
  }

  oldLen <- length(ms)
  nR1 <- nextn(oldLen, 2)
  if (oldLen != nR1) {
    ms <- c(ms, ms[oldLen:(2 * oldLen - nR1 + 1)])
  }

  len <- length(ms)

  wCoefs <- matrix(0,
                   ncol = length(scales), nrow = length(ms),
                   dimnames = list(1:length(ms), scales)
  )
  for (scale.i in scales) {
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax)) / (scale.i * dxval))
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    fft <- Re(fft(fft(ms) * Conj(fft(f)), inverse = TRUE)) / len
    wCoefs.i <- 1 / sqrt(scale.i) * fft
    ind <- len - floor(lenWave / 2)
    wCoefs[, as.character(scale.i)] <- c(wCoefs.i[(ind + 1):len], wCoefs.i[1:(ind)])
  }
  wCoefs[1:oldLen, , drop = FALSE]
}

doBinning <- function(spectra, split = "scan", tolerance = 0.002) {
  df_list <- split.data.frame(spectra, spectra[, split])

  nonEmpty <- sapply(df_list, nrow) != 0L # checking if the list is not empty
  samples <- rep.int(
    seq_along(df_list),
    sapply(df_list, nrow)
  )

  mass <- unname(unlist((lapply(
    df_list[nonEmpty],
    function(x) as.double(x$mz)
  )),
  recursive = FALSE, use.names = FALSE
  ))
  intensities <- unlist(lapply(
    df_list[nonEmpty],
    function(x) as.double(x$i)
  ),
  recursive = FALSE, use.names = FALSE
  )
  snr <- unlist(lapply(
    df_list[nonEmpty],
    function(x) x$SNR
  ),
  recursive = FALSE, use.names = FALSE
  )
  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]

  mass <- binning(
    mass = mass, intensities = intensities,
    samples = samples, tolerance = tolerance
  )

  s <- sort.int(mass, index.return = TRUE) # sort results into mass lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]
  lIdx <- split(seq_along(mass), samples)

  df_list[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses
    p <- NULL
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    p$snr <- snr[i]
    as.data.frame(p)
  }, p = df_list[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df_list
}

#' @title Bin spectra with similar mass ranges
#' @param peakList
#' @param fraction
#' @param npeaks
#' @param meanSNR
#' @param tolerance
#' @importFrom pbapply pblapply
#' @export
binSpectra <- function(peakList, fraction = 0, npeaks = 0,
                       meanSNR = 0, tolerance = 0.002) {
  pbapply::pblapply(peakList, function(spectra) {
    df <- do.call(rbind, doBinning(spectra, split = "scan", tolerance))
    df <- df[order(df$mz), ]
    bins <- split.data.frame(df, df$mz)
    df <- data.frame(
      mz = unique(df$mz), npeaks = sapply(bins, nrow),
      i = sapply(bins, function(x) sum(x$i)),
      SNR = sapply(bins, function(x) mean(x$snr))
    )
    rownames(df) <- 1:nrow(df)
    df[df$npeaks >= (fraction * length(unique(spectra$scan))) &
         df$SNR >= meanSNR & df$npeaks > npeaks, ]
  })
}

#' @title Bin cells
#' @param spectraList List of spectra
#' @param tolerance
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
binCells <- function(spectraList, sampleData = NULL,
                     tolerance = 0.002, filter = TRUE) {
  df <- do.call(rbind, lapply(1:length(spectraList), function(i) {
    df <- spectraList[[i]]
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$rt <- -1
    df$cell <- i
    df
  }))
  df$mz <- round(df$mz, 4)

  df_list <- doBinning(df, split = "cell", tolerance = tolerance)
  df <- data.frame(do.call(rbind, df_list),
                   cell = rep(names(df_list), sapply(df_list, nrow))
  )

  bins <- split.data.frame(df, df$mz)

  m <- matrix(NA, nrow = length(bins), ncol = length(df_list))
  snr <- matrix(NA, nrow = length(bins), ncol = length(df_list))

  colnames(m) <- as.integer(unique(df$cell))
  rownames(m) <- 1:nrow(m)
  dimnames(snr) <- dimnames(m)
  for (i in 1:length(bins)) {
    m[i, bins[[i]]$cell] <- bins[[i]]$intensity
    snr[i, bins[[i]]$cell] <- bins[[i]]$snr
  }
  SummarizedExperiment(
    colData = data.frame(sample = names(spectraList)),
    assays = list(Area = m, SNR = snr),
    rowData = data.frame(mz = as.double(names(bins)))
  ) %>% addSampleData(sampleData) %>% filterCells()
}
