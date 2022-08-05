file <- list.files(file.path(r"(F:\Myrthe\mzml)"), full.names = T)[2]
z <- mzR::openMSfile(file)
x <- mzR::peaks(z)
df <- mzR::header(z)

range <- abs(df$scanWindowUpperLimit - df$scanWindowLowerLimit)
scans <- range >= 0 & range <= Inf
polarity <- "+"
if (!is.null(polarity)) {
  polarity_filter <- grepl("FTMS - p NSI", df$filterString)
  if (polarity == "+") polarity_filter <- !polarity_filter
  scans <- scans & polarity_filter
}

spectrum <- x[scans]
spectrum <- x[[2]]

intensities <- sumR:::savgol(spectrum[, 2], 2)
intensities[intensities < 0] <- 0

## find local maxima
centroids <- which(diff(sign(diff(intensities))) == -2) + 1

mzInd <- spectrum[, 1]
noise <- sumR:::cwt_new(spectrum[, 2], scales = 1)
minNoiseLevel <- max(noise) * 0.0001
noise <- abs(noise)
nMz <- length(mzInd)
noiseWindow <- 0.01
winSize.noise <- as.integer(nMz * noiseWindow)


peakNoise <- sapply(centroids, function(k){
  start.k <- max(1, k - winSize.noise)
  end.k <- min(nMz, k + winSize.noise)
  noiseLevel.k <- quantile(noise[start.k:end.k], probs = 0.05, na.rm = T)
  if (noiseLevel.k < minNoiseLevel) {
    noiseLevel.k <- minNoiseLevel
  }
  noiseLevel.k
})

peakNoise <- calculateNoise(spectrum, centroids, noiseWindow)
spectrum <- as.data.frame(spectrum[centroids, ])
colnames(spectrum) <- c("mz", "i")
spectrum$Noise <- peakNoise
spectrum

