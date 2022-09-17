saverImputation <- function(df, cores = 1, normalize = TRUE, ...){
  df[is.na(df)] <- 0
  if (!normalize) normalize <- NULL
  if (requireNamespace("SAVER", quietly = TRUE)) {
    df <- suppressMessages(suppressWarnings(
      SAVER::saver(df, estimates.only = T, ncores = cores, size.factor = normalize
      )))
  } else {
    warning("Please install the package 'SAVER' before using this imputation")
  }
  df

}

magicImputation <- function(df, ...){
  df[is.na(df)] <- 0
  if (requireNamespace("Rmagic", quietly = TRUE)) {
    loadNamespace("Rmagic")
    if (Rmagic::pymagic_is_available()) {
      df <- t(Rmagic::magic(t(df))$result)
    }
  }
  df
}

rfImputation <- function(df, ...){
  if (requireNamespace("pmp", quietly = TRUE)) {
    loadNamespace("pmp")
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
#' @noRd
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
