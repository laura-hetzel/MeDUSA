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
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = phenotype)) +
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    theme(legend.position="top")}
  local.save_plot(paste("PCA",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
}

# *** Welch [T-Test] -----------------------------------------------------
#' MZLOG-OBJ T-Test
#'
#' Welch T-test, input should be an log2 of mz_obj
#'  - Requires: dplyr
#''
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param phenotype_a \cr
#'   DataFrame: Filtered metadata for phenotypeA
#' @param phenotype_b \cr
#'   DataFrame: Filtered metadata for phenotypeB
#' @param adjust_method \cr
#'   String: How to p.adjust?
#'   Boolean: False = no adjustment ; True = "fdr"
#'
#' Dependencies : pbapply, dplyr
#' @examples
#' To t.test values from "phenoA" over "phenoB"
#'
#' phenotype_a <- filter(metadata, phenotype == "phenoA" &
#'                       ionmode == ionmode_val & phase == phase_val &
#'                       filtered_out == "no")
#' phenotype_b <- filter(metadata, phenotype == "phenoB &
#'                       ionmode == ionmode_val  & phase == phase_val &
#'                       filtered_out == "no")
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
mzlog_analysis_welch <- function(input_mzlog_obj, phenotype_a, phenotype_b, adjust = 'fdr', cores = 2){
  phenotype_a <- phenotype_a
  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    out <- data.frame(p    = rep(Inf, nrow(input_mzlog_obj)),
                      p_05 = rep(FALSE, nrow(input_mzlog_obj)),
                      p_10 = rep(FALSE, nrow(input_mzlog_obj)),
                      p_15 = rep(FALSE, nrow(input_mzlog_obj)),
                      row.names = input_mzlog_obj$mz)

    out$p <- pbapply::pbapply(input_mzlog_obj,1 , cl = cl, function(i){
        t.test(dplyr::select(input_mzlog_obj[i,], phenotype_a$sample_name),
               dplyr::select(input_mzlog_obj[i,], phenotype_b$sample_name))$p.value
    })

    out$p_05 <- out$p < 0.05
    out$p_10 <- out$p < 0.10
    out$p_15 <- out$p < 0.15
    out <- tibble::rownames_to_column(out, "mz")
    out$mz <- as.numeric(out$mz)
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
  out$mz <- row.names(out)
  out
}


# *** Fold [change] -----------------------------------------------------
#' MZLOG-OBJ Fold Change
#'
#' Compare two phenotypes
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param phenotype_a \cr
#'   DataFrame: Filtered metadata for phenotypeA
#' @param phenotype_b \cr
#'   DataFrame: Filtered metadata for phenotypeB
#'
#' @examples
#' To fold values from "phenoA" over "phenoB"
#'
#' phenotype_a <- filter(metadata, phenotype == "phenoA" &
#'                       ionmode == ionmode_val & phase == phase_val &
#'                       filtered_out == "no")
#' phenotype_b <- filter(metadata, phenotype == "phenoB &
#'                       ionmode == ionmode_val  & phase == phase_val &
#'                       filtered_out == "no")
#' @param fold_math
#'   Method: How to combine Phenotype intensities
#'
#' @returns DataFrame with columns:
#'  - mz   : Float
#'  - fold : Float
#'
#' @export
mzlog_analysis_fold <- function(input_mzlog_obj, phenotype_a, phenotype_b, fold_math = "mean"){
  pheno_a <- input_mzlog_obj[,phenotype_a$sample_name]
  pheno_b <- input_mzlog_obj[,phenotype_b$sample_name]
  input_mzlog_obj$fold <- (apply(pheno_a, 1, fold_math) /
                           apply(pheno_b, 1, fold_math))
  return(dplyr::select(input_mzlog_obj, mz,fold))
}
