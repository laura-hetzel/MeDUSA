# *** PCA -----------------------------------------------------
#' MZLOG-OBJ PCA
#'
#' Not really sure
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
#' Welch T-test, input should be an log2 of mz_obj
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
    out <- tibble::rownames_to_column(out, "mz")
    out$mz <- as.numeric(df_l$mz)
    return(out)
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
  #out$mz <- row.names(out)
  out
}


# *** Fold [change] -----------------------------------------------------
#' MZLOG-OBJ Fold Change
#'
#' Compare two phenotypes
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
