#' @title Imputation of missing values using SAVER
#' @inherit imputation description details
#' @returns A data.frame with all values filled with the same dimensions
#' as the given argument in `data`.
#' @param data Data.frame of the assay to be imputated
#' @param normalize Boolean value, should data be normalized before imputation?
#' @inherit imputation examples
saverImputation <- function(data, normalize = TRUE){
  data[is.na(data)] <- 0
  if (!normalize) normalize <- NULL
  if (requireNamespace("SAVER", quietly = TRUE)) {
    data <- suppressMessages(suppressWarnings(
      quiet(SAVER::saver(data, estimates.only = T, ncores = 1,
                         size.factor = normalize)
      )))
  } else {
    warning("Please install the package 'SAVER' before using this imputation")
  }
  data

}

#' @title Imputation of missing values using Rmagic
#' @inherit imputation description details
#' @returns A data.frame with all values filled with the same dimensions
#' as the given argument in `data`.
#' @param data Data.frame of the assay to be imputated
#' @inherit imputation examples
magicImputation <- function(data){
  data[is.na(data)] <- 0
  if (requireNamespace("Rmagic", quietly = TRUE)) {
    loadNamespace("Rmagic")
    if (Rmagic::pymagic_is_available()) {
      data <- t(quiet(Rmagic::magic(t(data))$result))
    }
  } else {
    warning("Please install the package 'Rmagic' before using this imputation")
  }
  data
}

#' @title Imputation of missing values using models
#' @inherit imputation description details
#' @returns A data.frame with all values filled with the same dimensions
#' as the given argument in `data`.
#' @param data Data.frame of the assay to be imputated
#' @param model Character, name of the model to be used for missing value
#' imputation. Requires the `pmp` packages and must be one of the
#' methods mentioned in `mv_imputation`.
#' @inherit imputation examples
modelImputation <- function(data, model = "rf"){
  if (requireNamespace("pmp", quietly = TRUE)) {
    loadNamespace("pmp")
    data <- quiet(pmp::mv_imputation(data, method = model))
  } else {
    warning("Please install the package 'pmp' before using this imputation")
  }
  data
}

#' @title Imputate missing values using random noise
#' @inherit imputation description details
#' @returns A data.frame with all values filled with the same dimensions
#' as the given argument in `data`.
#' @param data Data.frame of the assay to be imputated
#' @param noise numerical value for the maximum random noise level. Defaults
#' to 100.
#' @importFrom stats runif
#' @inherit imputation examples
noiseImputation <- function(data, noise = 100) {
  idx <- is.na(data) | data == 0
  data[idx] <- runif(sum(idx), min = 1, max = noise)
  return(data)
}

#' @title Impute missing values of an assay
#' @description This function aids in imputation of peaks that were not found
#' with [extractPeaks]. Currently, four methods are supported of imputation
#' are supported through the `method` parameter. These are based on
#' random (noise) or model imputation ([modelImputation], [magicImputation],
#' and [saverImputation]). The former uses regular models like KNN or Random
#' Forest models to impute values. The latter two are based on single-cell RNA
#' methods that are more complex and take correlation between cells into
#' account.
#' @details Imputation of missing values is a vital part of post-processing
#' as it is needed for other post-processing methods, including statistics and
#' modelling.
#' @returns A SummarizedExperiment object with an added slot in assays given
#' with the `saveAssay` parameter.
#' @param exp SummarizedExperiment object obtained after alignment, containing
#' at least the slot given in `useAssay`.
#' @param method Character vector of method to use for imputation. Can be one
#' of `"noise"`, `"model"`, `"magic"`, `"saver"`. Defaults to `"noise"`.
#' @param useAssay Character name of the assay to be imputed. Defaults to
#' `"Area"`
#' @param saveAssay Character name of the to be saved imputed assay. Defaults
#' to `"Imputed"`
#' @param setDefault Boolean value, should the imputed assay be set at index 1,
#' making it the default wher using `assay(exp)`? Defaults to `TRUE`
#' @inheritDotParams noiseImputation -data
#' @inheritDotParams modelImputation -data
#' @inheritDotParams saverImputation -data
#' @importFrom SummarizedExperiment assay<-
#' @export
#' @examples
#' # Loading negative polarity dataset
#' data("sumRnegative")
#'
#' # Noise Imputation
#' exp <- imputation(sumRnegative, method = "noise", noise = 250)
#'
#' # Model Imputation
#' exp <- imputation(sumRnegative, method = "model", model = "sv")
#'
#' # Rmagic Imputation
#' # exp <- imputation(sumRnegative, method = "magic")
#'
#' # SAVER Imputation
#' exp <- imputation(sumRnegative, method = "saver", normalize = TRUE)
imputation <- function(exp, method = "noise", useAssay = "Area",
                       saveAssay = "Imputed", setDefault = TRUE, ...) {

  if (!validateExperiment(exp)) return(NULL)
  exp <- filterCells(exp)
  df <- as.data.frame(assay(exp, useAssay))

  new_df <- switch(method[1],
                   "saver" = saverImputation(df, ...),
                   "magic" = magicImputation(df),
                   "model" = modelImputation(df, ...),
                   "noise" = noiseImputation(df, ...)
  )
  dimnames(new_df) <- dimnames(assay(exp))
  assay(exp, saveAssay) <- new_df
  if (setDefault) exp <- setDefaultAssay(exp, saveAssay)
  metadata(exp)$Imputation <- method[1]
  exp

}
