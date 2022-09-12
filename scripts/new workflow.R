library(sumR)

folder <- r"(C:\Users\Pascal Maas\Downloads\dominique)"
mzmls <- rawToMzml(folder, file.path(folder, "mzML"))


peaks <- extractPeaks(
  mzmls,
  polarity = "+",
  cores = 8
)


peaks <- filterScansByBiomarker(peaks, mass = 368.42)
spectra <- spectraAlignment(peaks, method = "binning", ppm = 5, nPeaks = 3)
spectra
exp <- cellAlignment(spectra,  method = "binning", ppm = 5,  nCells = 5)
spectra

peak <- peaks[["sample21_220907143140"]][as.integer(unlist(spectra[22821, "peakIdx"])), ]
max(rowData(exp)$ppm)

peaks[[1]]
peaks[[2]]
lapply(peaks, function(a){
  a[a$mz >= 992.446 & a$mz <= 992.454, ]
})

