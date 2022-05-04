pickSpectra <- function(x, SNR=0){
  ms <- as.data.frame(x)
  colnames(ms) <- c("mz", "i")
  #    return(performPeakPicking(ms, SNR))
  wCoefs <- cbind("0" = ms$i, cwt_new(ms$i))
  localMax <- MassSpecWavelet::getLocalMaximumCWT(wCoefs)
  ridgeList <- MassSpecWavelet::getRidge(localMax)
  majorPeakInfo <- identifyPeaks(ms$i, ridgeList, wCoefs,
                                 SNR.Th = SNR, nearbyPeak = TRUE)

  df <- cbind(ms[as.vector(majorPeakInfo$allPeakIndex), ],
              Noise = majorPeakInfo$peakNoise[names(majorPeakInfo$allPeakIndex)],
              Value = majorPeakInfo$peakValue[names(majorPeakInfo$allPeakIndex)],
              SNR = majorPeakInfo$peakSNR[names(majorPeakInfo$allPeakIndex)])
  colnames(df) <- make.names(colnames(df), unique = T)
  df[df$SNR >= SNR, ]
}

peakPicking <- function(files, SNR = 0){
  result <- pbapply::pblapply(files, function(f){
    z <- mzR::openMSfile(f)
    x <- mzR::peaks(z)
    l <- lapply(seq_len(length(x)), function(i){
      tryCatch({
        suppressWarnings(pickSpectra(x[[i]], SNR))
      }, error = function(err) NULL)
    })

    non_nulls <- !vapply(l, is.null, logical(1))
    l <- l[non_nulls]
    df <- data.frame(scan = rep(which(non_nulls), sapply(l, nrow)), do.call(rbind, l))
    dplyr::distinct(df[df$Value > 0, ])
  })
  samps <- tools::file_path_sans_ext(basename(f))

  names(result) <- samps
  result
}

identifyPeaks <- function (ms, ridgeList, wCoefs, scales = as.numeric(colnames(wCoefs)),
                           SNR.Th = 3, peakScaleRange = 5, ridgeLength = 32, nearbyPeak = FALSE,
                           nearbyWinSize = 100, winSize.noise = 500, SNR.method = "quantile",
                           minNoiseLevel = 0.001) {

  if (ridgeLength > max(scales))
    ridgeLength <- max(scales)
  peakScaleRange <- scales[scales >= peakScaleRange]

  minNoiseLevel <- max(wCoefs) * minNoiseLevel

  ridgeLen <- sapply(ridgeList, length)
  ridgeName <- names(ridgeList)
  ridgeInfo <- matrix(as.numeric(unlist(strsplit(ridgeName,
                                                 "_"))), nrow = 2)
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
    }
    else {
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
    ms.int <- ms[start.k:end.k]
    noiseLevel.k <- switch(SNR.method, quantile = quantile(noise[start.k:end.k],
                                                           probs = 0.95), sd = sd(noise[start.k:end.k]), mad = mad(noise[start.k:end.k],
                                                                                                                   center = 0), data.mean = mean(ms.int), data.mean.quant = mean(ms.int[ms.int <
                                                                                                                                                                                          quantile(ms.int, probs = 0.95)]))
    if (noiseLevel.k < minNoiseLevel)
      noiseLevel.k <- minNoiseLevel
    peakSNR <- c(peakSNR, peakValue[k]/noiseLevel.k)
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
  selInd3 <- !(mzInd %in% c(1:(nearbyWinSize/2), (nrow(wCoefs) -
                                                    (nearbyWinSize/2) + 1):nrow(wCoefs)))
  selInd <- (selInd1 & selInd2 & selInd3)
  names(peakSNR) <- names(peakScale) <- names(peakNoise) <- names(peakCenterInd) <- names(peakValue) <- names(mzInd) <- ridgeName
  return(list(peakIndex = mzInd[selInd], peakValue = peakValue,
              peakCenterIndex = peakCenterInd, peakSNR = peakSNR,
              peakNoise = peakNoise,
              peakScale = peakScale,
              potentialPeakIndex = mzInd[selInd1 & selInd3], allPeakIndex = mzInd))
}


cwt_new <- function(ms, scales=NULL, max_scale = 32){
  psi_xval <- seq(-8, 8, length = 1024)
  psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2/2)

  psi_xval <- psi_xval - psi_xval[1]
  dxval <- psi_xval[2]
  xmax <- psi_xval[length(psi_xval)]

  if (is.null(scales)){
    scales <- seq(1, ifelse(floor(length(ms) / 16) < max_scale,
                            floor(length(ms) / 16), max_scale), 2)
  }

  oldLen <- length(ms)
  nR1 <- nextn(oldLen, 2)
  if (oldLen != nR1) {
    ms <- c(ms, ms[oldLen:(2 *  oldLen - nR1 + 1)])
  }

  len <- length(ms)

  wCoefs <- matrix(0, ncol = length(scales), nrow = length(ms),
                   dimnames = list(1:length(ms), scales))
  for (scale.i in scales) {
    f <- rep(0, len)
    j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
    lenWave <- length(j)
    f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
    fft <- Re(fft(fft(ms) * Conj(fft(f)), inverse = TRUE)) / len
    wCoefs.i <- 1/sqrt(scale.i) * fft
    ind <- len - floor(lenWave/2)
    wCoefs[, as.character(scale.i)] <- c(wCoefs.i[(ind + 1):len], wCoefs.i[1:(ind)])
  }
  wCoefs[1:oldLen, , drop = FALSE]
}

groupPeaks <- function(l, fraction = 0.1, meanSNR = 0.5){
  result <- pbapply::pblapply(l, function(df){
    df$sample <- df$scan
    df$rt <- -1
    x <- suppressMessages(xcms::do_groupChromPeaks_density(peaks = df,
                                          sampleGroups = as.factor(rep(1, nrow(df))),
                                          minFraction = 0,
                                          binSize = 0.025))

    x$i <- vapply(1:nrow(x), function(i) sum(df[unlist(x[i, "peakidx"]), "i"]), numeric(1))

    x$meanSNR <- vapply(1:nrow(x), function(i) mean(df[unlist(x[i, "peakidx"]), "SNR"]), numeric(1))
    x <- x[, -which(colnames(x) == "peakidx")]
    x[x$npeaks >= fraction * length(unique(df$scan)) & x$meanSNR >= meanSNR, ]
  })
  names(result) <- names(l)
  result
}


groupCells <- function(res){
  df <- do.call(rbind, lapply(1:length(res), function(i){
    df <- res[[i]]
    if (nrow(df) == 0) return(NULL)
    df$rt <- -1
    df$mz <- df$mzmed
    df$sample <- i
    df
  }))
  x <- xcms::do_groupChromPeaks_density(peaks = df,
                                        sampleGroups = as.factor(names(res)),
                                        minFraction = 0,
                                        binSize = 0.25)

  m <- matrix(nrow = nrow(x), ncol = length(res))
  for (i in 1:nrow(x)) {
    ids <- unlist(x[i, "peakidx"])
    m[i, df$sample[ids]] <- df$i[ids]
  }
  dimnames(m) <- list(rownames(x), names(res))
  x <- x[, -which(colnames(x) == "peakidx")]

  SummarizedExperiment(assays = m, rowData = x)
}


library(pbapply)
library(SummarizedExperiment)
library(xcms)
library(magrittr)

dir <- file.path(r"(C:\Users\pmaas\Documents\GitHub\sum-r\sumR)")
f <- list.files(dir, full.names = T, pattern = ".mzML")
f
#data <- readMSData(f, mode = "onDisk")
#data <- smooth(data)
#data <- pickPeaks(data)
#data <- data[fData(data)$scanWindowUpperLimit - fData(data)$scanWindowLowerLimit <= 200, ]
#output <- sprintf("%s_centroided.mzML", tools::file_path_sans_ext(basename(f)))
#writeMSData(data, output)



exp <- peakPicking(f) %>%
  groupPeaks(fraction = 0.1, meanSNR = 0.5) %>%
  groupCells()


exp
assay(exp)

rowData(exp)



library(Rmagic)
library(SAVER)

df <- as.data.frame(assay(exp[,1:10]))
df[is.na(df)] <- 0

MAGIC_data <- magic(df)

sav <- saver(df)

head(MAGIC_data$result)[,1:5]
head(sav$estimate)[,1:5]
head(df)[,1:5]

pca <- prcomp(df)
p1 <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point() +
  ggtitle("raw pca")

pca <- prcomp(MAGIC_data$result)
p2 <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point() +
  ggtitle("imputed pca MAGIC")

pca <- prcomp(sav$estimate)
p3 <- ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point() +
  ggtitle("imputed pca SAVER")

cowplot::plot_grid(plotlist = list(p1, p2, p3), ncol = 1)
