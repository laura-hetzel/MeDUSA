# library(sumR)
# filenames <- list.files("C:/Users/Pascal Maas/Downloads/", pattern = ".*.mzML", full.names = T)
# file_df <- as.data.frame(do.call(rbind, lapply(filenames, function(mzml){
#   # mzml <- "C:/Users/Pascal Maas/Downloads/Blank-EC-LOW-1.mzML"
#   combine_spectra_centroid(read_mzml(mzml))
# })))
#
# library(sumR)
# dir <- system.file("extdata", package = "sumR")
# files <- list.files(dir, pattern = ".*negative.*.mzML$", full.names = T)
#
# data <- dims_pipeline(read_mzml(files))
#
# df <- feature_processing(data)
# df <- isotope_tagging(df)
