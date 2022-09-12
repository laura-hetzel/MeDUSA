
# BG filter done by Amar

#' @title Blank filter to remove background
#' @description It is used to conduct background removal, it filters a data frame of m/z values and intensities
#' @param dataframe a data.frame of m/z values and intensities including blank samples
#' @param limit intensity ratio limit, for example within 150 percent of the blank intensity = 1.5
#' @param blank_regex regex of blank columns, for example "blank"
#' @param sample_regex regex of sample columns for example
#' @usage bgFilter(dataframe, limit = 2.5, blank_regex = "blank", sample_regex = "^d|^ec|qc")
#' @export
bgFilter <- function(dataframe, limit = 2.5, blank_regex = "blank", sample_regex = "^d|^ec|qc") {
  ## Extract blank column names
  blank_cols <- grep(blank_regex, names(dataframe), ignore.case = TRUE)

  ## Extract sample column names
  sample_cols <- grep(sample_regex, names(dataframe), ignore.case = TRUE)

  ## Create new column with intensity means of blank samples
  dataframe$background <- rowMeans(dataframe[, blank_cols], na.rm = TRUE)

  ## Compare sample_cols with the mean, replace them by 0 if they are below
  ## thresh from what I understand up until now: Sweep creates the Intensity
  ## Ratio as a summary statistic, which sweeps across the sample columns
  ## The function across those columns is basically a division by the
  ## corresponding Background value.
  ## I'm gettin errors related to: Dimensions of matrix need to be equal to the
  ## function input data
  ## https://stackoverflow.com/questions/3444889/how-to-use-the-sweep-function

  ## Make a subset with intensity columns
  sample_df <- dataframe[, sample_cols]

  ## Filter the intensity columns by conditioning to smaller or equal to limit
  ## Returns 0 if condition is met
  sample_df <- as.data.frame(lapply(sample_df, function(i) {
    ifelse(i / dataframe$background <= limit, 0, i)
  }))

  # merys : needs to be checked if return is right
  return(sample_df)

  # dataframe[sample_cols][sweep(dataframe[sample_cols], 1,
  # dataframe$background, `/`) <= limit] <- 0
  # result <- dataframe[, -blank_cols]
  # return(result)
  # ahmeds notes: Need to review this later, maybe replace swap with map_dfc function
}

#' MD filter
#'
#' @description filters data based on Mass Defect
#'
#' @importFrom dplyr mutate select
#'
#' @param dataframe is the dataframe on which KM calculations need to be
#' conducted
#' @param mz_col is the mass to charge column within the specified dataframe
#' @param a is the coefficient of the linear equation used for data filtering
#' @param b is the addition value of the linear equation used for data
#' filtering, both a and b can be calculated using the linear_equation()
#' function included in this package.
#'
#' @usage MD_filter(dataframe, mz_col, a = 0.00112, b = 0.01953)
#' @export
#' @importFrom dplyr filter %>%
#' @import ggplot2
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
      spectra[spectra$scan %in% scans[2:(length(scans) - 1)], ]
    })
}

#' @title filter peaks
#' @export
peakFilter <- function(fileList, intensity = 0, peakValue = -Inf, snr = 0){
  lapply(fileList, function(spectra){
    intensity <- spectra$i >= intensity
    peakValue <- spectra$i - spectra$Noise > peakValue
    snr <- spectra$i / spectra$Noise > snr # snr <- (spectrum$i - spectrum$Noise) / spectrum$Noise > snr
    spectra[intensity & peakValue & snr, ]
  })
}

#' @title filter spectra
#' @export
spectraFilter <- function(spectra, npeaks = 0, intensity = 0, SNR = 0){
  res <- lapply(spectra, function(x) x[x$SNR >= SNR & x$i >= intensity & x$npeaks >= npeaks, ])
  attributes(res) <- attributes(spectra)
  res
}


setDefaultAssay <- function(exp, default){
  if (!validateExperiment(exp)) return(NULL)
  n <- assayNames(exp)
  to_replace <- which(n == default)
  n[to_replace] <- n[1]
  n[1] <- default
  assays(exp) <- assays(exp)[n]
  exp
}

saverImputation <- function(df, cores = 1, normalize = TRUE, ...){
  df[is.na(df)] <- 0
  if (!normalize) normalize <- NULL
  if (!"SAVER" %in% installed.packages()) {
    warning("Please install the package 'SAVER' before using this imputation")
  } else {
    df <- suppressMessages(suppressWarnings(
      SAVER::saver(df, estimates.only = T, ncores = cores, size.factor = normalize
      )))
  }
  df

}

magicImputation <- function(df, ...){
  df[is.na(df)] <- 0
  if ("Rmagic" %in% installed.packages()) {
    if (Rmagic::pymagic_is_available()) {
      df <- t(Rmagic::magic(t(df))$result)
    }
  }

  df
}

rfImputation <- function(df, ...){
  if ("pmp" %in% installed.packages()) {
    df <- pmp::mv_imputation(df, method = "rf")
  }
  df
}

#' @title data imputation
#' @description This function apply data imputation from 1 to selected noise level
#' @param data dataframe
#' @param noise numerical value for the noise level
#' @param seed global seed for reproducible results
#' @importFrom stats runif
noiseImputation <- function(data, noise = 100, seed = 42, ...) {
  set.seed(seed)
  idx <- is.na(data) | data == 0
  data[idx] <- runif(sum(idx), min = 1, max = noise)
  return(data)
}

#' @title Impute missing values
#' @importFrom SummarizedExperiment assay<-
#' @export
imputation <- function(exp, method = "noise", useAssay = "Area",
                       saveAssay = "Imputed", setDefault = TRUE, ...) {

  if (!validateExperiment(exp)) return(NULL)
  exp <- filterCells(exp)
  df <- as.data.frame(assay(exp, useAssay))

  new_df <- switch(method,
    "saver" = saverImputation(df, ...),
    "magic" = magicImputation(df, ...),
    "rf" = rfImputation(df, ...),
    "noise" = noiseImputation(df, ...)
  )
  dimnames(new_df) <- dimnames(assay(exp))
  assay(exp, saveAssay) <- new_df
  if (setDefault) exp <- setDefaultAssay(exp, saveAssay)
  exp

}

#' @title Fragment Filter
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

