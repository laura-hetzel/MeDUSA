# *** PCA -----------------------------------------------------
#' MZ-OBJ PCA
#'
#' Not really sure
#'  - Requires: ggplot2, tibble
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param metadata \cr
#'   DataFrame?: MZ_metadata
#' @param bad_samples \cr
#'   List?     : c("bad1", "bad2")
#' @param high_noise \cr
#'   Float     : HighBoundry for "Noise"
#'
#' @export
mzlog_pca <- function(input_mzlog_obj,metadata, bad_samples ) {
  t_mz_obj <- scale(t(input_mzlog_obj))
  t_mz_obj <- t_mz_obj[!row.names(t_mz_obj) %in% bad_samples ]
  t_mz_obj <- merge(x = t_mz_obj,
        y = subset(tibble::column_to_rownames(metadata, "sample_name"), select = "phenotype"),
        by = 0, all.x = TRUE)
  rownames(t_mz_obj) <- paste(t_mz_obj$Row.names, t_mz_obj$phenotype, sep = "&")

  t_mz_obj <- prcomp(subset(t_mz_obj, select = -c(Row.names, phenotype)))
  summary(t_mz_obj)

  var_explained <- t_mz_obj$sdev^2/sum(t_mz_obj$sdev^2)

  pca_neg$x %>%
    as.data.frame %>%
    rownames_to_column("sample_phenotype") %>%
    separate(sample_phenotype, c("sample", "phenotype"), "&") %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(color = phenotype)) +
    labs(x=paste0("PC1: ",round(var_explained_n[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained_n[2]*100,1),"%")) +
    theme(legend.position="top")
  local.save_plot(paste("PCA",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
}

# *** T-Test -----------------------------------------------------
#' MZ-OBJ T-Test
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
#'
#' @examples
#' To remove values from "samples" that are lower than "blanks" * "threshold"
#'
#' phenotype_a <- filter(metadata, phenotype == "phenoA" &
#'                       ionmode == ionmode_val & phase == phase_val &
#'                       filtered_out == "no")
#' phenotype_b <- filter(metadata, phenotype == "phenoB &
#'                       ionmode == ionmode_val  & phase == phase_val &
#'                       filtered_out == "no")
#'
#' @returns DataFrame with columns:
#'  - rowname: Float: MZ
#'  - p      : Float: p-value
#'  - p_05   : Boolean: is P < 0.05
#'  - p_10   : Boolean: is P < 0.10
#'  - p_15   : Boolean: is P < 0.15
#'
#' @export
mzlog_welch <- function(input_mzlog_obj, phenotype_a, phenotype_b, cores = 4){
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, varlist = ls(environment(mzlog_welch)),
                              envir = environment(mzlog_welch))
  } else {
    cl <- NULL
  }
  tryCatch({
    out <- data.frame(p    = rep(Inf, nrow(input_mzlog_obj)),
                      p_05 = rep(FALSE, nrow(input_mzlog_obj)),
                      p_10 = rep(FALSE, nrow(input_mzlog_obj)),
                      p_15 = rep(FALSE, nrow(input_mzlog_obj)),
                      row.names = row.names(input_mzlog_obj))

    out$p <- pbapply::pblapply(1:nrow(input_mzlog_obj), cl = cl, function(i){
        t.test(dplyr::select(input_mzlog_obj[i,], phenotype_a$sample_name),
               dplyr::select(input_mzlog_obj[i,], phenotype_b$sample_name))$p.value
    })

    out$p_05 <- out$p < 0.05
    out$p_10 <- out$p < 0.10
    out$p_15 <- out$p < 0.15
    return(out)
  },
  finally={
    if (cores > 1) {
      parallel::stopCluster(cl)
      showConnections()
    }
  })
}

# *** RandomForest -----------------------------------------------------
#' MZ-OBJ RandomForest
#'
#' Not really sure
#'  - Requires: ggplot2, tibble
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Ob
#' @param metadata \cr
#'   DataFrame?: MZ_metadata
#' @param bad_samples \cr
#'   List?     : c("bad1", "bad2")
#' @param high_noise \cr
#'   Float     : HighBoundry for "Noise"
#'
#' @export
mzlog_random_forest <- function(input_mzlog_obj){}
