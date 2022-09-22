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

#' @title Plot chromatogram of a peak for 1 or more samples
#' @description This function plots the chromatogram of an aligned peak. This
#' can be used to inspect the peak alignment.
#' @param exp SummarizedExperiment object obtained after alignment with
#' [cellAlignment].
#' @param peaks List of centroided peaks obtained from [extractPeaks]
#' @param sampleIdx Vector of 1 or more sample indices or names to plot the
#' chromatogram for.
#' @param peakId Peak index or peak name to plot the chromatogram for.
#' @importFrom stats aggregate.data.frame
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw .data
#' @export
chromatogramPlot <- function(exp, peaks, sampleIdx = 1, peakId = 1){
  centroided <- peaks[sampleIdx]

  peak <- rowData(exp[peakId, ])
  df <- do.call(rbind, lapply(sampleIdx, function(i){
    df <- centroided[[i]]
    df <- df[df$rt >= peak$rtmin & df$rt <= peak$rtmax &
                    df$mz >= peak$mzmin & df$mz <= peak$mzmax, ]
    df <- aggregate(df, by = list(rt = df$rt), FUN = "sum")
    df$sample <- i
  }))

  ggplot(as.data.frame(df), aes(x = .data$rt, y = .data$i,
                                color = .data$sample)) +
    geom_point() +
    geom_line() +
    theme_bw()
}




#' @title spectraShiftPlot
#' @param spectraList
#' @param file File number of the files that were used.
#' @importFrom ggplot2 ggplot geom_segment aes unit xlab arrow
#' @export
spectraShiftPlot <- function(spectra, cell = 1){
  if (is.numeric(cell)){
    cell <- unique(spectra$sample)[cell]
  }

  if (!cell %in% spectra$sample) {
    stop(sprintf("Cell '%s' not found in spectra", cell))
  }

  spectra$mzmax <- spectra$mzmax - spectra$mz
  spectra$mzmin <- spectra$mzmin - spectra$mz
  ggplot(as.data.frame(spectra[spectra$sample == cell, ])) +
    geom_segment(aes(x = mzmin, xend = mzmax, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    geom_segment(aes(x = mzmax, xend = mzmin, y = mz, yend = mz),
                 arrow = arrow(length = unit(2, "mm"))) +
    xlab("mz Shift")
}


#' @title cellShiftPlot
#' @param exp
#' @importFrom ggplot2 ggplot geom_segment aes unit xlab arrow
#' @export
cellShiftPlot <- function(exp) {
  if (!validateExperiment(exp)) return(NULL)

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

