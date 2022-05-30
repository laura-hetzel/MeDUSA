#' @title Plot spectrum peaks
#' @export
plotSpectrumPeaks <- function(peaks, file = 1, scan = 1, SNR = 0){
  peaks <- peakFilter(peaks, SNR)
  file <- attr(peaks, "Files")[file]
  df <- as.data.frame(sumR:::formatScans(file,
                                         attr(peaks, "massWindow"),
                                         attr(peaks, "polarity"),
                                         attr(peaks, "combineSpectra"))[[scan]])

  colnames(df) <- c("mz", "i")
  wCoefs <- cbind("0" = df$i, sumR:::cwt_new(df$i))

  stacked <- stack(wCoefs)
  stacked$mz <- df[stacked$row, "mz"]
  stacked <- as.data.frame(stacked)

  low <- stacked[stacked$col == colnames(wCoefs)[2], ]

  base <- tools::file_path_sans_ext(basename(file))
  result <- peaks[[base]]
  result <- result[result$scan == scan, ]
  if (!is.null(df) & nrow(result) > 0) {
    ggplot(df) +
      geom_point(data = result, aes(x = mz, y = i, color = "Peaks")) +
      geom_path(data = low, aes(x = mz, y = value, color = "Wavelet")) +
      geom_line(data = result, aes(x = mz, y = Noise, color = "Noise"), linetype = "dashed") +
      geom_segment(aes(x = mz, y = 0, yend = i, xend = mz))
  }
}

#' @title noisePlot
#' @export
noisePlot <- function(peaks, file = 1){
  result <- dplyr::distinct(peaks[[file]][,c("scan", "Noise")])
  ggplot(result, aes(x = scan, y = Noise)) +
    geom_line() +
    ggplot2::scale_y_log10()
}


#' @title plotCellPeaks
#' @export
plotCellPeaks <- function(peaks, file = 1){
  ggplot(peaks[[file]], aes(x = mz, y = scan, size = SNR)) +
    geom_density_2d_filled() +
    geom_point(shape = 18) +
    theme_minimal()
}

#' @title plotSpectraShift
#' @export
plotSpectraShift <- function(spectraList, whichSpectrum = 1){
  ggplot(as.data.frame(spectraList[[whichSpectrum]])) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz Shift")
}


#' @title plotCellShift
#' @export
plotCellShift <- function(exp) {
  ggplot(as.data.frame(rowData(exp))) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz shift")
}
