#' @title Pick peaks in a spectrum
#' @importFrom MassSpecWavelet getLocalMaximumCWT getRidge
pickSpectrum <- function(x, SNR = 0, scales = NULL, maxScale = 32, noiseWindow = 0.1,
                        solventFilter = NULL) {
  ms <- as.data.frame(x)
  colnames(ms) <- c("mz", "i")
  if (!is.null(solventFilter)){
    solventFilter <- solventFilter[solventFilter$mzmin >= min(ms$mz) & solventFilter$mzmax <= max(ms$mz), ]

    ms$i <- sapply(1:nrow(ms), function(i){
      check <- any(ms[i, "mz"] >= solventFilter$mzmin & ms[i, "mz"] <= solventFilter$mzmax)
      if (check){
        return(0)
      }
      return(ms[i, "i"])
    })
  }
  winSize.noise <- abs(max(ms$mz) - min(ms$mz)) * noiseWindow
  wCoefs <- cbind("0" = ms$i, cwt_new(ms$i, scales, max_scale = maxScale))
  localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs, minWinSize = 0)
  ridgeList <- MassSpecWavelet::getRidge(localMax, minWinSize = 0)
  majorPeakInfo <- identifyPeaks(ms$i, ridgeList, wCoefs, winSize.noise = winSize.noise,
                                 SNR.Th = SNR, nearbyPeak = F, peakScaleRange = 0
  )

  df <- cbind(ms[as.vector(majorPeakInfo$allPeakIndex), ],
              Noise = majorPeakInfo$peakNoise[names(majorPeakInfo$allPeakIndex)],
              Value = majorPeakInfo$peakValue[names(majorPeakInfo$allPeakIndex)],
              SNR = majorPeakInfo$peakSNR[names(majorPeakInfo$allPeakIndex)]
  )
  colnames(df) <- make.names(colnames(df), unique = T)
  df
}

#' @title Perform peak picking of a single file
#' @importFrom pbapply pblapply
#' @export
#'
peakPickFile <- function(file = NULL, scales = c(1,4,9), noiseWindow = 0.1,
                         spectra = NULL, ...){

  if (is.null(spectra)) spectra <- prepareFile(file, ...)

  logging::loginfo("Peak picking %s", basename(file))
  l <- pbapply::pblapply(spectra, function(spectrum) {
    tryCatch(suppressWarnings(pickSpectrum(spectrum, scales = scales,
                                          noiseWindow = noiseWindow
                                          )),
    error = function(err) NULL)
  })

  non_nulls <- !vapply(l, is.null, logical(1))
  l <- l[non_nulls]
  df <- data.frame(scan = rep(which(non_nulls), sapply(l, nrow)), do.call(rbind, l))

  if (nrow(df) == 0) return(NULL)
  df <- df[order(df$mz), ]
  rownames(df) <- 1:nrow(df)
  df
}

#' @title Perform peak picking of a single file
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @export
peakPickFiles <- function(files, ...){
  result <- lapply(setNames(files, basename(files)), function(file) peakPickFile(file, ...) )
  result <- result[!vapply(result, is.null, logical(1))]

  attributes(result) <- list(...)
  result
}

#' @title Perform peak picking
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores makeCluster clusterExport stopCluster
#' @importFrom mzR openMSfile peaks header
#' @importFrom dplyr distinct
#' @importFrom tools file_path_sans_ext
#' @export
peakPicking <- function(files, maxScale = 32, noiseWindow = 0.1, cores = 1,
                        solventFilter = NULL, scales = NULL) {
  cl <- NULL
  if (cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = names(sys.frame()))
  }

  result <- lapply(1:length(files), function(i) {
    logging::loginfo("Peak picking %s", names(files)[i])
    l <- pbapply::pblapply(fileList[[i]], cl = cl, function(spectrum) {

      tryCatch(suppressWarnings(pickSpectra(spectrum, scales = scales,
                                            maxScale = maxScale,
                                            noiseWindow = noiseWindow,
                                            solventFilter = solventFilter)),
               error = function(err) NULL)
    })
    non_nulls <- !vapply(l, is.null, logical(1))
    l <- l[non_nulls]
    df <- data.frame(scan = rep(which(non_nulls), sapply(l, nrow)), do.call(rbind, l))

    if (nrow(df) == 0) return(NULL)
    df <- df[order(df$mz), ]
    rownames(df) <- 1:nrow(df)
    df
  })
  if (!is.null(cl)) stopCluster(cl)

  result <- result[!vapply(result, is.null, logical(1))]
  attributes(result) <- attributes(fileList)
  result
}

identifyPeaks <- function(ms, ridgeList, wCoefs, scales = as.numeric(colnames(wCoefs)),
                          SNR.Th = 3, peakScaleRange = 5, ridgeLength = 32, nearbyPeak = FALSE,
                          nearbyWinSize = 100, winSize.noise = 50, SNR.method = "quantile",
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
  scales <- as.numeric(colnames(wCoefs))
  noise <- abs(wCoefs[, which.min(scales[scales > 0])]) #"1"])
  peakSNR <- NULL
  peakNoise <- NULL
  nMz <- nrow(wCoefs)
  for (k in 1:length(ridgeList)) {
    ind.k <- mzInd[k]
    start.k <- max(1, ind.k - winSize.noise)
    end.k <- min(nMz, ind.k + winSize.noise)

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
  val <- 8
  psi_xval <- seq(-val, val, length = 1024)
  psi <- (2 / sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2 / 2)

  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax <- psi_xval[length(psi_xval)]

  if (is.null(scales)) {
    scales <- seq(1, ifelse(floor(length(ms) / 16) < max_scale,
                            floor(length(ms) / 16), max_scale), 2)
  }

  oldLen <- length(ms)
  nR1 <- nextn(oldLen, 2)
  if (oldLen != nR1) {
    ms <- c(ms, ms[oldLen:(2 * oldLen - nR1 + 1)])
  }

  len <- length(ms)

  wCoefs <- matrix(0, ncol = length(scales), nrow = length(ms),
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

