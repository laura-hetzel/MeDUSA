#' @title Plot spectrum peaks
#' @export
spectrumPlot <- function(peaks, file = 1, scan = 1){
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
cellPlot <- function(peaks, file = 1){
  ggplot(peaks[[file]], aes(x = mz, y = scan, size = SNR)) +
    geom_density_2d_filled() +
    geom_point(shape = 18) +
    theme_minimal()
}

#' @title spectraShiftPlot
#' @export
spectraShiftPlot <- function(spectraList, file = 1){
  ggplot(as.data.frame(spectraList[[file]])) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz Shift")
}


#' @title cellShiftPlot
#' @export
cellShiftPlot <- function(exp) {
  df <- as.data.frame(rowData(exp))
  df$mzmin <- df$mzmin - df$mz
  df$mzmax <- df$mzmax - df$mz
  ggplot(df) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz shift")
}

#' @title Plot feature coverage when adding cells
#' @export
featureCovPlot <- function(cells, assay = 1, by = 5, seed = 42){
  set.seed(seed)
  cellNames <- sample(colnames(cells))
  is <- seq(by, ncol(cells), by)
  feats <- vapply(is, function(i){
    m <- assay(cells[, cellNames[1:i]], assay)
    nrow(m) - sum(rowSums(is.na(m)) == ncol(m))
  }, double(1))
  df <- data.frame(Cells = is, Features = feats)
  ggplot(df, aes(x = Cells, y = Features)) +
    geom_line()
}

#' @title Plot features per cell
#' @export
featureCellPlot <- function(cells, assay = "Area"){
  df <- as.data.frame(stack(colSums(!is.na(assay(cells, assay)))))
  df <- df[order(df$values, decreasing = T), ]
  df$Cell <- 1:nrow(df)
  colnames(df)[1] <- "Features"
  ggplot(df, aes(y = Features, x = Cell)) + geom_bar(stat = "identity")
}
