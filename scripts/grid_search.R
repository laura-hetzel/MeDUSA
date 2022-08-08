dir <- file.path(r"(C:\Users\pmaas\Downloads\Myrthe\Mythe)")
files <- list.files(path = dir, pattern = ".mzML$", full.names = TRUE)

library(sumR)
counts <- lapply(setNames(files, basename(files)), function(file){
  sum(sapply(mzR::peaks(mzR::openMSfile(file)), function(df){
    any(df[,1] > 372.21 & df[,1] < 372.25)
  }))
})

counts_profile <- do.call(rbind, counts)
counts <- as.data.frame(counts_profile)


fileList <- prepareScans(files[c(8, 10, 12)], polarity = "+", centroid = FALSE)

counts_centroided <- lapply(fileList, function(spectra){
  sum(sapply(spectra, function(df){
    any(df[,1] > 372.21 & df[,1] < 372.25)
  }))
})
counts$centroided <- unlist(counts_centroided)
counts
colnames(counts) <- c('Profile', "Positive Polarity", "Centroided")



#### Blank substraction

fileList <- peakPicking(fileList, cores = 3, maxScale = 8)

spectraList <- spectraBinning(fileList)
bins <- cellBinning(spectraList)



grid <- expand.grid(
  scales = list(NULL, c(1, 4, 9), seq(1:4)),
  maxScale = c(4, 8, 16, 32),
  combineSpectra = c(FALSE, TRUE)
)

windows <- grid$noiseWindow
scals <- grid$scales
comb <- grid$combineSpectra
i <- 1
peaks <- peakPicking(files,
                     noiseWindow = windows[i],
                     scales = unlist(scals[i]))


