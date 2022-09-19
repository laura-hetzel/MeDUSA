#' @title Substract Blanks from SummarizedExperiment
#' @description
#' @details
#' @returns
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
#' @examples
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

#' @title Filter cells without measurements (only containing NA values)
#' @description Not all cells will contain the same peaks after alignment. For
#' peaks that are not found in cells, an `NA` value is used instead. This
#' function checks if the entire cell does not only contain NAs and remove
#' them if so.
#' @details During post-processing, several filters may be applied to remove
#' either compounds or peaks to reduce false positive hits. In some cases,
#' cells with few peaks may not contain any peaks at all. This function
#' helps in detecting and removing such cells.
#' @returns
#' @param exp SummarizedExperiment object obtained after alignment
#' @param assay The assay to be checked. Defaults to the assay at index 1
#' @export
#' @examples
filterCells <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  exp[, colSums(is.na(assay(exp, assay))) != nrow(exp)]
}

#' @title Filter scans based on Biomarkers like oil or solvent
#' @description This function can be used to extract data obtained from the
#' "slug" method that is currently in development. It will extract all scans
#' after the first mass with intensity > threshold is found and stop when that
#' same mass is found again.
#' @details The "slug" method results in a total ion chromatogram consisting
#' of three parts. The first and last part contain ionization solvent at
#' higher concentrations than the second part, which contains the cell
#' measurement. By extracting this middle part from the total ion chromatogram,
#' the only remaining part is the cell with measurements of interest.
#' @param fileList List of centroided peaks obtained from `extractPeaks()`
#' @param mass Biomarker mass to use as identification in the solvent
#' @param intensity Intensity value that is expected in the solvent
#' @param ppm Mass error for the biomarker. Defaults to 20
#' @export
#' @examples
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

#' @title Identify and remove fragments from aligned peaks
#' @description This function finds correlation based on intensities given
#' in the assay. Peaks that have a higher correlation than the value
#' given in `corr` will be considered fragments and are removed.
#' @details Fragments occur in mass spectrometry quite often. Since these
#' fragments originate from the same compound / metabolite, the resulting
#' intensities are expected to be very similar. By using either `spearman` or
#' `pearson` correlation, the correlation between peaks can be calculated
#' and be removed.
#'
#' If a lot of values were originally missing, it is recommended to use
#' spearman correlation instead of pearson correlation as the latter tends to
#' give higher correlation values. Spearman correlation uses ranked based
#' correlation and has shown to give more accurate correlation values when
#' many values are imputed.
#' @returns
#' @param exp SummarizedExperiment object obtained after alignment
#' @param assay Name or index of the assay to use. Defaults to the first
#' assay (index 1)
#' @param method One of either "spearman" or "pearson". It is used as method
#' in the stats `cor` function.
#' @param corr Minimum correlation for peaks to be considered fragments.
#' Defaults to 0.95
#' @importFrom stats cor
#' @export
#' @examples
fragmentFilter <- function(exp, assay = 1, method = "spearman", corr = 0.95){
  if (!validateExperiment(exp)) return(NULL)
  corrs <- cor(t(assay(exp, assay)), method = method)
  max_corr <- vapply(1:nrow(corrs), function(i) max(corrs[i, -i]), double(1))
  exp[max_corr < corr, ]
}

#' @title Remove grouped peaks (features) per phenotype
#' @description
#' @details
#' @returns
#' @param exp SummarizedExperiment obtained after alignment
#' @param assay Assay to be used, defaults to the assay at index 1
#' @param nCells Minimum number of cells a feature should be found in. Defaults
#' to 0, meaning all features are kept.
#' @param fCells Similar to `nCells` but for fractions instead. Can be used
#' to remove features that are found in less than half of the cells. Defaults
#' to 0, meaning all features are kept.
#' @export
#' @examples
featureFilter <- function(exp, assay = 1, nCells = 0, fCells = 0){
  if (!validateExperiment(exp)) return(NULL)

  to_take <- do.call(intersect, lapply(unique(metadata(exp)$phenotype),
                                       function(split){
    sub <- exp[, exp$phenotype == split]
    nFeats <- as.double(apply(assay(sub, assay), 1, function(x) sum(!is.na(x))))
    pFeats <- nFeats / ncol(sub)
    which(nFeats >= nCells & pFeats >= fCells)
  }))
  exp[to_take, ]
}

#' @title Keep the most variable features based on Levene Test
#' @description This function removes features that have similar features
#' across the given phenotypes and calculated by the function `leveneTest`.
#' @details When interested in which peaks are significantly different between
#' groups in the given phenotypes, the leveneTest can be used. Here, only
#' features are kept where the variance between the groups is significantly
#' different. This greatly reduces the dataset with often only a few compounds
#' remaining.
#' @returns A SummarizedExperiment with only compounds that have different
#' variances across the phenotype given.
#' @param exp SummarizedExperiment object after alignment and leveneTest
#' @export
#' @examples
keepVariableFeatures <- function(exp){
  if (!validateExperiment(exp)) return(NULL)
  if (!"leveneTest" %in% colnames(rowData(exp))) {
    message("'leveneTest' not found in rowData. Please run 'leveneTest' first.")
    return(exp)
  }
  exp[which(rowData(exp)$leveneTest$unequal_variance), ]
}

