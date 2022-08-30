#' @title Plot spectrum peaks
#' @param peaks List of datafres with peak picked data
#' @param file File number of the files that were used.
#' @param scan
#' @importFrom ggplot2 .data
#' @export
spectrumPlot <- function(peaks, file = 1, scan = 1){
  df <- peaks[[file]]
  df <- df[df$scan == unique(df$scan)[scan], ]
  ggplot(df, aes(x = .data$mz, y = .data$i)) +
    geom_segment(aes(x = .data$mz, y = 0, yend = .data$i, xend = .data$mz))
}

#' @title noisePlot over all spectra
#' @param peaks List of datafres with peak picked data
#' @param file File number of the files that were used.
#' @export
noisePlot <- function(peaks, file = 1){
  result <- dplyr::distinct(peaks[[file]][,c("scan", "Noise")])
  ggplot(result, aes(x = scan, y = Noise)) +
    geom_line() +
    ggplot2::scale_y_log10()
}


#' @title plotCellPeaks
#' @param peaks List of datafres with peak picked data
#' @param file File number of the files that were used.
#' @export
cellPlot <- function(peaks, file = 1){
  ggplot(peaks[[file]], aes(x = mz, y = scan, size = SNR)) +
    geom_density_2d_filled() +
    geom_point(shape = 18) +
    theme_minimal()
}

#' @title spectraShiftPlot
#' @param spectraList
#' @param file File number of the files that were used.
#' @export
spectraShiftPlot <- function(spectra, cell = 1){
  if (is.numeric(cell)){
    cell <- unique(spectra$sample)[cell]
  }

  if (!cell %in% spectra$sample) stop(sprintf("Cell '%s' not found in spectra", cell))

  spectra$mzmax <- spectra$mzmax - spectra$mz
  spectra$mzmin <- spectra$mzmin - spectra$mz
  ggplot(as.data.frame(spectra[spectra$sample == cell, ])) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz), arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz Shift")
}


#' @title cellShiftPlot
#' @param exp
#' @export
cellShiftPlot <- function(exp) {
  if (!validateExperiment(exp, checkColData = F)) return(NULL)

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
#' @param cells
#' @param assay
#' @param by
#' @param seed
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
  loess_fit <- loess(Features ~ Cells, df)
  df2 <- as.data.frame(predict(loess_fit))
  df2$Cells <- df$Cells
  colnames(df2) <- c("Features", "Cells")

  ggplot(df, aes(x = Cells, y = Features)) +
    geom_line(aes(colour = "Features")) +
    geom_line(data = df2, aes(colour = "Trend")) +
    labs(colour = "Line Type")
}

#' @title Plot features per cell
#' @param cells
#' @param assay
#' @export
featureCellPlot <- function(cells, assay = "Area"){
  df <- as.data.frame(stack(colSums(!is.na(assay(cells, assay)))))
  df <- df[order(df$values, decreasing = T), ]
  df$Cell <- 1:nrow(df)
  colnames(df)[1] <- "Features"
  ggplot(df, aes(y = Features, x = Cell)) + geom_bar(stat = "identity")
}

#' @title Plot feature spectra
#' @param cells
#' @param feature
#' @param file File number of the files that were used.
#' @export
plotRawFeature <- function(cells, feature, file){
  xmin <- rowData(cells)$mzmin[feature]
  xmax <- rowData(cells)$mzmax[feature]

  x <- formatScans(file, 200, "+", F)
  x2 <- do.call(rbind, lapply(x, function(df){
    df[df[,1] >= xmin & df[,1] <= xmax, ]
  }))
  if (nrow(x2) == 0) return(NULL)
  plot3D::scatter3D(x = x2[, 1], y = 1:nrow(x2), z = x2[, 2], theta = 45, phi = 10,
            bty = "g",  type = "h", ylab = "Scan",
            xlab = "mz", zlab = "i",
            ticktype = "detailed", pch = 19, cex = .5)
}
