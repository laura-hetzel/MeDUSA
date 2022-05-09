doBinning <- function(spectra, split = "scan", tolerance = 0.002){
  df_list <- split.data.frame(spectra, spectra[,split])

  nonEmpty <- sapply(df_list, nrow) != 0L #checking if the list is not empty
  samples <- rep.int(seq_along(df_list),
                     sapply(df_list, nrow))

  mass <- unname(unlist((lapply(df_list[nonEmpty],
                                function(x) as.double(x$mz))),
                        recursive = FALSE, use.names = FALSE))
  intensities <- unlist(lapply(df_list[nonEmpty],
                               function(x) as.double(x$i)),
                        recursive = FALSE, use.names = FALSE)
  snr <- unlist(lapply(df_list[nonEmpty],
                       function(x) x$SNR),
                recursive = FALSE, use.names = FALSE)
  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]

  mass <- binning(mass = mass, intensities = intensities,
                  samples = samples, tolerance = tolerance)

  s <- sort.int(mass, index.return = TRUE)# sort results into mass lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  snr <- snr[s$ix]
  lIdx <- split(seq_along(mass), samples)

  df_list[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses
    p = NULL
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
binSpectra <- function(peakList, fraction = 0, npeaks = 0,
                       meanSNR = 0, tolerance = 0.002){
  pbapply::pblapply(peakList, function(spectra){
    df <- do.call(rbind, doBinning(spectra, split = "scan", tolerance))
    df <- df[order(df$mz), ]
    bins <- split.data.frame(df, df$mz)
    df <- data.frame(mz = unique(df$mz), npeaks = sapply(bins, nrow),
                     i = sapply(bins, function(x) sum(x$i)),
                     SNR = sapply(bins, function(x) mean(x$snr)))
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
binCells <- function(spectraList, tolerance = 0.002){
  df <- do.call(rbind, lapply(1:length(spectraList), function(i){
    df <- spectraList[[i]]
    if (nrow(df) == 0) return(NULL)
    df$rt <- -1
    df$cell <- i
    df
  }))
  df$mz <- round(df$mz, 4)

  df_list <- doBinning(df, split = "cell", tolerance = tolerance)
  df <- data.frame(do.call(rbind, df_list),
                   cell = rep(names(df_list), sapply(df_list, nrow)))

  bins <- split.data.frame(df, df$mz)

  m <- matrix(NA, nrow = length(bins), ncol = length(df_list))
  snr <- matrix(NA, nrow = length(bins), ncol = length(df_list))

  colnames(m) <- as.integer(unique(df$cell))
  rownames(m) <- 1:nrow(m)
  dimnames(snr) <- dimnames(m)
  for (i in 1:length(bins)) {
    m[i, bins[[i]]$cell] <- bins[[i]]$intensity
    snr[i, bins[[i]]$cell] <- bins[[i]]$snr
  }
  SummarizedExperiment(assays = list(Area = m, SNR = snr),
                       rowData = data.frame(mz = as.double(names(bins))))
}

#' Filter cells without measurements
#' @export
filterCells <- function(exp){
  exp[, colSums(is.na(assay(exp))) != nrow(exp)]
}



#' @importFrom SAVER saver
#' @importFrom SummarizedExperiment assay<-
#' @export
imputation <- function(exp, normalize = TRUE, useAssay = "Area",
                       saveAssay = "Imputed", cores = 1){
  df <- as.data.frame(assay(exp, useAssay))
  df[is.na(df)] <- 0
  assay(exp, saveAssay) <- saver(df, estimates.only = T, ncores = cores,
                                 size.factor = as.integer(normalize))
  exp
}
