# *** PCA -----------------------------------------------------
#' MZLOG-OBJ PCA
#'
#' A PCA plot is a relatively simple and well received method to compare samples
#' and determine if they can be grouped by phenotype. The mzlog_analysis_pca 
#' function scales the log transformed data, removes any blacklisted m/z from
#' consideration, and incorporates the metadata to create a PCA plot that is 
#' colored by phenotype for easy interpretation.
#' 
#'  - Requires: ggplot2, tibble
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param metadata \cr
#'   DataFrame: MZ_metadata
#' @param sample_blacklist \cr
#'   List?     : c("bad1", "bad2")
#'
#' @export
mzlog_analysis_pca <- function(input_mzlog_obj,metadata, sample_blacklist = c() ) {
  metadata <- local.meta_polarity_fixer(input_mzlog_obj,metadata)
  rownames(input_mzlog_obj) <- input_mzlog_obj$mz
  t_mz_obj <- scale(t(dplyr::select(input_mzlog_obj,-mz)))
  t_mz_obj <- t_mz_obj[!row.names(t_mz_obj) %in% sample_blacklist ,]
  t_mz_obj <- merge(x = t_mz_obj,
        y = subset(tibble::column_to_rownames(metadata, "sample_name"), select = "phenotype"),
        by = 0, all.x = TRUE)
  rownames(t_mz_obj) <- paste(t_mz_obj$Row.names, t_mz_obj$phenotype, sep = "&")

  t_mz_obj <- prcomp(subset(t_mz_obj, select = -c(Row.names, phenotype)))
  summary(t_mz_obj)

  var_explained <- t_mz_obj$sdev^2/sum(t_mz_obj$sdev^2)

  t_mz_obj$x %>%
    as.data.frame %>%
    rownames_to_column("sample_phenotype") %>%
    separate(sample_phenotype, c("sample", "phenotype"), "&") %>%

    ggpubr::ggscatter(x = "PC1",y="PC2",
                      title = paste("PCA: ", local.mz_polarity_guesser(input_mzlog_obj)),
                      color ="phenotype",
                      size  = 0.5,
                      ellipse = TRUE,
                      mean.point = TRUE,
                      mean.point.size = 2,
                      xlab = paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                      ylab = paste0("PC2: ",round(var_explained[2]*100,1),"%"))

  local.save_plot(paste("PCA",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
}

# *** Welch [T-Test] -----------------------------------------------------
#' MZLOG-OBJ T-Test
#'
#' To determine if there is a statistical difference between phenotypes, the 
#' Welch T-rest is often used. The user is expected to isolate two phenotypes in
#' the data (example provided for code to isolate) to pass into the 
#' mzlog_analysis_welch function, which uses Welch t-test function and an FDR
#' correction to output the p-vale for each m/z.
#' 
#'  - Requires: dplyr
#''
#' @param phenoA_mz_obj \cr
#'   DataFrame: mz_obj of phenoA samples
#' @param phenoB_mz_obj \cr
#'   DataFrame: mz_obj of phenoB samples
#' @param adjust_method \cr
#'   String: How to p.adjust?
#'   Boolean: False = no adjustment ; True = "fdr"
#'
#' Dependencies : pbapply, dplyr
#' @examples
#' To t.test values from "phenoA" over "phenoB"
#'
#' phenoA_mz_obj  <- sumR::mztools_filter(input_mzObj,metadata,"phenoA")
#' phenoB_mz_obj  <- sumR::mztools_filter(input_mzObj,metadata,"phenoB")
#'
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#'
#' @returns DataFrame with columns:
#'  - mz      : Float: MZ
#'  - p       : Float: p-value
#'  - p_05    : Boolean: is P < 0.05
#'  - p_10    : Boolean: is P < 0.10
#'  - p_15    : Boolean: is P < 0.15
#'  - adjusted: p.adjust($p, adjust_method) *optional, but default
#'
#' @export
mzlog_analysis_welch <- function(phenoA_mz_obj, phenoB_mz_obj, adjust = 'fdr', cores = 2){
  df_l <- local.ensure_mz(phenoA_mz_obj,phenoB_mz_obj, "sumR::mzlog_analysis_welch")

  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    out <- data.frame(p    = rep(Inf, nrow(df_l$mz)),
                      p_05 = rep(FALSE, nrow(df_l$mz)),
                      p_10 = rep(FALSE, nrow(df_l$mz)),
                      p_15 = rep(FALSE, nrow(df_l$mz)))

    out$p <- pbapply::pbapply(df_l$mz, 1 , cl = cl, function(i){
      t.test(df_l$df_a[df_l$mz==i,],
             df_l$df_b[df_l$mz==i,])$p.value
    })

    out$p_05 <- out$p < 0.05
    out$p_10 <- out$p < 0.10
    out$p_15 <- out$p < 0.15
    return(cbind(df_l$mz, out))
  },
  finally={
    local.kill_threads(cl)
  })

  if(adjust != F){
    if (adjust){
      adjust = 'fdr'
    }
    out$adjusted <- as.data.frame(p.adjust(out$p), method = adjust_method)
  }
  out
}


# *** Fold [change] -----------------------------------------------------
#' MZLOG-OBJ Fold Change
#'
#' Fold change is one method of comparing two phenotypes, dividing one phenotype
#' by the other. The user is expected to isolate two phenotypes in the data 
#' (example provided for code to isolate) to pass into the mzlog_analysis_fold 
#' function, which will calculate the mean intensity (default, user defined) of 
#' each phenotype per m/z. The mean of one phenotype is then divided by the 
#' mean of the other phenotype, yielding one unitless fold change value for 
#' each m/z. This allows the user to determine which m/z values are highly 
#' up/down regulated or mildly up/down regulated.
#'
#' @param phenoA_mz_obj \cr
#'   DataFrame: mz_obj of phenoA samples
#' @param phenoB_mz_obj \cr
#'   DataFrame: mz_obj of phenoB samples
#' @param fold_math
#'   Method: How to combine Phenotype intensities
#'
#' @examples
#' To fold values from "phenoA" over "phenoB"
#'
#' phenoA_mz_obj  <- sumR::mztools_filter(input_mzObj,metadata,"phenoA")
#' phenoB_mz_obj  <- sumR::mztools_filter(input_mzObj,metadata,"phenoB")
#'
#' @returns DataFrame with columns:
#'  - mz   : Float
#'  - fold : Float
#'#'

#' @export
mzlog_analysis_fold <- function(phenoA_mz_obj, phenoB_mz_obj, fold_math = "mean"){
  df_l <- local.ensure_mz(phenoA_mz_obj,phenoB_mz_obj, "sumR::mzlog_analysis_fold")
  out <- df_l$mz
  out$fold <- (apply(df_l$df_a, 1, fold_math) /
               apply(df_l$df_b, 1, fold_math))
  return(out)
}
