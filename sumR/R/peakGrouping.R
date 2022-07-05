doBinning <- function(spectra, split = "scan", tolerance = 0.002) {
  df_list <- split.data.frame(spectra, spectra[, split])

  nonEmpty <- sapply(df_list, nrow) != 0L # checking if the list is not empty
  samples <- rep.int(
    seq_along(df_list),
    sapply(df_list, nrow)
  )

  mass <- unname(unlist((lapply(
    df_list[nonEmpty],
    function(x) as.double(x$mz)
  )),
  recursive = FALSE, use.names = FALSE
  ))
  intensities <- unlist(lapply(
    df_list[nonEmpty],
    function(x) as.double(x$i)
  ),
  recursive = FALSE, use.names = FALSE
  )
  snr <- unlist(lapply(
    df_list[nonEmpty],
    function(x) x$SNR
  ),
  recursive = FALSE, use.names = FALSE
  )
  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]

  mass <- binning(
    mass = mass, intensities = intensities,
    samples = samples, tolerance = tolerance
  )

  s <- sort.int(mass, index.return = TRUE) # sort results into mass lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]
  lIdx <- split(seq_along(mass), samples)

  df_list[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses
    p <- NULL
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    p$snr <- snr[i]
    as.data.frame(p)
  }, p = df_list[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  df_list
}

#' @title Bin spectra with similar mass ranges
#' @param peakList
#' @param fraction
#' @param npeaks
#' @param meanSNR
#' @param tolerance
#' @importFrom pbapply pblapply
#' @export
spectraBinning <- function(peakList, fraction = 0, npeaks = 0, method = "sum",
                           meanSNR = 0, tolerance = 0.002) {
  pbapply::pblapply(peakList, function(spectra) {
    if (nrow(spectra) == 0) return(NULL)
    df_list <- doBinning(spectra, split = "scan",
                         tolerance = tolerance)

    df <- do.call(rbind, df_list)
    if (nrow(df) == 0) return(NULL)

    df <- df[order(df$mz), ]
    df$oldmz <- spectra$mz[order(spectra$mz)]

    df$mzdiff <- df$mz - df$oldmz
    df <- df[order(rownames(df)), ]
    bins <- split.data.frame(df, df$mz)
    df <- data.frame(
      mz = unique(df$mz), npeaks = sapply(bins, nrow),
      i = sapply(bins, function(x) switch(method, "sum" = sum(x$i), "mean" = mean(x$i))),
      SNR = sapply(bins, function(x) mean(x$snr)),
      mzmin = sapply(bins, function(x) min(x$mzdiff)),
      mzmax = sapply(bins, function(x) max(x$mzdiff))
    )

    rownames(df) <- 1:nrow(df)
    df[df$npeaks >= (fraction * length(unique(spectra$scan))) &
         df$SNR >= meanSNR & df$npeaks > npeaks, ]
  })
}

#' @title Bin cells
#' @param spectraList List of spectra
#' @param tolerance
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
cellBinning <- function(spectraList, cellData = NULL, phenotype = NULL,
                        tolerance = 0.002, filter = TRUE) {
  df <- do.call(rbind, lapply(1:length(spectraList), function(i) {
    df <- spectraList[[i]]
    if (is.null(df)) return(NULL)
    if (nrow(df) == 0) return(NULL)

    df$rt <- -1
    df$cell <- names(spectraList)[i]
    df
  }))
  df$mz <- round(df$mz, 4)

  df_list <- doBinning(df, split = "cell", tolerance = tolerance)


  df2 <- do.call(rbind, df_list)
  df2 <- df2[order(df2$mz), ]
  df2$oldmz <- df$mz[order(df$mz)]
  df2$mzdiff <- df2$mz - df2$oldmz
  df2 <- df2[order(rownames(df2)), ]


  df <- data.frame(df2, cell = rep(names(df_list), sapply(df_list, nrow)))
  bins <- split.data.frame(df, df$mz)

  m <- matrix(NA, nrow = length(bins), ncol = length(df_list))
  rownames(m) <- 1:nrow(m)
  colnames(m) <- unique(df$cell)

  snr <- matrix(NA, nrow = length(bins), ncol = length(df_list))
  dimnames(snr) <- dimnames(m)

  for (i in 1:length(bins)) {
    m[i, bins[[i]]$cell] <- bins[[i]]$intensity
    snr[i, bins[[i]]$cell] <- bins[[i]]$snr
  }

  missing <- names(spectraList)[!names(spectraList) %in% colnames(m)]
  m <- cbind(m, matrix(NA, nrow = nrow(m), ncol = length(missing),dimnames = list(rownames(m), missing)))
  snr <- cbind(snr, matrix(NA, nrow = nrow(snr), ncol = length(missing), dimnames = list(rownames(snr), missing)))

  mzs <- as.double(names(bins))
  se <- SummarizedExperiment(
    assays = list(Area = m, SNR = snr),
    rowData = DataFrame(mz = mzs,
                        mzmin = mzs + sapply(bins, function(x) min(x$mzdiff)),
                        mzmax = mzs + sapply(bins, function(x) max(x$mzdiff))
    ),
    metadata = list(phenotype = phenotype)
  )
  if (!is.null(cellData)) colData(se) <- DataFrame(cellData[colnames(m), ])
  se
}

