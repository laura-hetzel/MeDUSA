#' @title Substract Blanks from SummarizedExperiment
#' @param exp SummarizedExperiment obtained after alignment
#' @param foldChange Multiplier used for the median blanks. If the Area of
#' a sample is lower than this threshold, it counts towards the `nSamples`
#' threshold. (Defaults to 5).
#' @param sampleThresh Number of samples are allowed to be lower than the
#' median blank area times the number given for `foldChange`. Defaults to `Inf`,
#' meaning that no compounds will be removed unless changed.
#' @param removeBlanks Should samples of the type `BLANK` be removed after the
#' filter is applied? Defaults to `TRUE`.
#' @export
blankSubstraction <- function(exp, foldChange = 5, nSamples = Inf,
                              removeBlanks = TRUE){
  if (!validateExperiment(exp)) return(exp)
  blanks <- exp[, toupper(exp$Type) == 'BLANK']
  samps <- exp[, toupper(exp$Type) != 'BLANK']

  threshold <- rowMedians(assay(blanks, "Area"), na.rm = TRUE) * foldChange
  exp <- exp[rowSums(assay(samps, "Area") - threshold <= 0, na.rm = TRUE) <= nSamples, ]

  if (removeBlanks) exp <- exp[, toupper(exp$Type) != "BLANK"]

  filterCells(exp)
}

#' @title MD filter
#' @description filters data based on Mass Defect
#' @param dataframe is the dataframe on which KM calculations need to be
#' conducted
#' @param mz_col is the mass to charge column within the specified dataframe
#' @param a is the coefficient of the linear equation used for data filtering
#' @param b is the addition value of the linear equation used for data
#' filtering, both a and b can be calculated using the linear_equation()
#' function included in this package.
#' @importFrom dplyr filter %>% mutate select
#' @import ggplot2
#' @export
MD_filter <- function(dataframe, mz_col, a = 0.00112, b = 0.01953) {
  ## Solves problem of customizable lin. eq as well: In case HMDB updates data
  ## In-function MD calculation
  dataframe$MZ <- mz_col
  ## Either floor() or trunc() can be used for this part.
  MZR <- trunc(mz_col, digits = 0)
  dataframe$MD <- MZ - MZR
  dataframe$MD.limit <- b + a * mz_col


  # dataframe <- dataframe %>%
  #   dplyr::mutate(MD, MZ, MD.limit) %>%
  #   dplyr::select(MD, MZ, MD.limit)

  ## This approach uses the linear equation to create an additional column and
  ## based on this column, we will be filtering out data.

  ## Below are in total 3 different plots: 1 of the starting dataframe, the
  ## second of the filtered dataframe the third is basically added to the first,
  ## in order to highlight the "to be removed" datapoints as a means of tagging
  ## them prior to removal.


  ## Notice how this is the exact opposite from the "filtered" dataframe below.
  ## This one filters everything ABOVE the limit, as everything below will be
  ## kept and we want those datapoints which will be removed, highlighted.
  ## Not those which will remain.
  highlight_df <- dataframe %>% filter(MD >= MD.limit)

  ## DO the same for inlcusion list of HMDB, We need to make: 1, Linear eq. 2,
  ## inclusion list plot 3, metaboshiny. 4, R package.

  MD_plot <- ggplot(data = dataframe, aes(x = MZ, y = MD)) +
    geom_point() +
    geom_point(data = highlight_df, aes(x = MZ, y = MD), color = "red") + # I added this one, so the data which will be removed will be highlighted in red.
    ggtitle(paste("Unfiltered MD data")) # deparse(substitute(dataframe)) if you want to add the dataframe name, acc to St. Ovflw
  # stat_smooth(method="lm", se=FALSE)-> For linear line through the plot, but may not be necessary to show


  # Creating a filtered dataframe:
  filtered <- dataframe %>%
    filter(MD <= MD.limit) # As I understood: Basically all are coordinates. The maxima equation basically gives coordinates
  # for the m/z values (x and MD = y). If it exceeds the equivalent coordinate of Y (which is MD) for the linear equation, it will be filtered.

  MD_plot_2 <- ggplot(data = filtered, aes(x = MZ, y = MD)) + # Filtered is basically the second dataframe, #which subsets datapoints with an Y value (which is the MD), below the linear equation MD...
    geom_point() +
    ggtitle(paste("Filtered MD data"))
  # stat_smooth(method="lm", se=FALSE) -> For linear line through the plot, but may not be necessary to show
  N_Removed_datapoints <- nrow(dataframe) - nrow(filtered) # To determine the number of peaks removed
  print(paste("Number of peaks removed:", N_Removed_datapoints))
  MD_PLOTS <- list(
    preMD_df = dataframe,
    postMD_df = filtered,
    preMD_plot = MD_plot,
    postMD_plot = MD_plot_2
  )
  return(MD_PLOTS)
}

#' Filter cells without measurements
#' @export
filterCells <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  exp[, colSums(is.na(assay(exp, assay))) != nrow(exp)]
}

#' @title Filter scans based on Biomarkers like oil or solvent
#' @param fileList List of centroided peaks obtained from `extractPeaks()`
#' @param mass Biomarker mass to use as identification in the solvent
#' @param intensity Intensity value that is expected in the solvent
#' @param ppm Mass error for the biomarker. Defaults to 20
#' @export
filterScansByBiomarker <- function(fileList, mass, intensity = 1e5, ppm = 20){
    minMass <- mass - mass * 1e-6 * ppm
    maxMass <- mass + mass * 1e-6 * ppm
    lapply(fileList, function(spectra){
      solvent_idx <- which(spectra$mz > minMass & spectra$mz < maxMass &
                             spectra$i > intensity)

      gap <- which.max(diff(solvent_idx))
      start <- solvent_idx[gap]
      end <- solvent_idx[gap + 1]

      spectra <- spectra[start:end, ]
      scans <- unique(spectra$scan)
      spectra <- spectra[spectra$scan %in% scans[2:(length(scans) - 1)], ]
      rownames(spectra) <- 1:nrow(spectra)
      spectra
    })
}

#' @title filter spectra
#' @export
spectraFilter <- function(spectra, npeaks = 0, intensity = 0, SNR = 0){
  res <- lapply(spectra, function(x) x[x$SNR >= SNR & x$i >= intensity & x$npeaks >= npeaks, ])
  attributes(res) <- attributes(spectra)
  res
}


#' @title Fragment Filter
#' @importFrom stats cor
#' @export
fragmentFilter <- function(exp, assay = 1, method = "spearman", corr = 0.95){
  if (!validateExperiment(exp)) return(NULL)
  corrs <- cor(t(assay(exp, assay)), method = method)
  max_corr <- vapply(1:nrow(corrs), function(i) max(corrs[i, -i]), double(1))
  exp[max_corr < corr, ]
}

#' @title Feature filter
#' @export
featureFilter <- function(exp, assay = "Area", nCells = 0, pCells = 0){
  if (!validateExperiment(exp)) return(NULL)
  to_take <- do.call(intersect, lapply(unique(exp$phenotype), function(split){
    sub <- exp[, exp$phenotype == split]
    nFeats <- as.double(apply(assay(sub, assay), 1, function(x) sum(!is.na(x))))
    pFeats <- nFeats / ncol(sub)
    which(nFeats >= nCells & pFeats >= pCells)
  }))
  exp[to_take, ]
}

#' @title Keep the most variable features
#' @param exp
#' @param assay
#' @param top
#' @export
keepVariableFeatures <- function(exp){
  if (!validateExperiment(exp)) return(NULL)
  if (!"leveneTest" %in% colnames(rowData(exp))) {
    message("'leveneTest' not found in rowData. Please run 'leveneTest' first.")
    return(exp)
  }
  exp[which(rowData(exp)$leveneTest$unequal_variance), ]
}

