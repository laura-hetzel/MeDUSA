# *** PCA -----------------------------------------------------
#' MZLOG-OBJ PCA
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
#' @param adjust \cr
#'   Boolean: Do a p.adjust?
#' @param adjust_method \cr
#'   String: How to p.adjust?
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
#'  - rowname : Float: MZ
#'  - p       : Float: p-value
#'  - p_05    : Boolean: is P < 0.05
#'  - p_10    : Boolean: is P < 0.10
#'  - p_15    : Boolean: is P < 0.15
#'  - adjusted: p.adjust($p, adjust_method) *optional, but default
#'
#' @export
mzlog_welch <- function(input_mzlog_obj, phenotype_a, phenotype_b, adjust = T, adjust_method = "fdr", cores = 4){

  cl <- local.export_thread_env(cores, deparse(sys.calls()[[sys.nframe()]]))
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
    if (cores > 1 || !is.null(cl)) {
      parallel::stopCluster(cl)
      showConnections()
    }
  })

  if(adjust){
    out$adjusted <- as.data.frame(p.adjust(out$p), method = adjust_method)
  }

}


volcano_fold_change <- function(df, oxy_numerator = "A"){
  #df <- df
  uniq_mz <- unique(df$mz)
  df_out <- data.frame(mz = uniq_mz,
                       fold = numeric(length(uniq_mz)))
  if ( nrow(filter(df, intensity == 0)) > 0 ) {
    stop(paste("ERROR: ", deparse(substitute(df)), "has loads of zeros (can't divide)"))
  }
  #cl <- parallel::makeCluster(4)
  #parallel::clusterExport(cl, varlist = ls(environment(fold_change)),
  #                        envir = environment(fold_change))
  df_out$fold <- pbapply::pbapply(df_out, 1, function(row){
    tmp_mz <- filter(df, mz == as.numeric(row[1]))
    tmp_mz$type <- gsub(".*_0[0-9]+_([A-Z]).*_(.*)", "\\1-\\2", tmp_mz$sample)
    tmp_squash <- aggregate(intensity ~ type , tmp_mz, mean)
    tmp_squash <- cbind(tmp_squash, str_match(tmp_squash$type, "(?<letter>[A-Z])-(?<cell>.*)"))

    fold <- ((filter(tmp_squash, type == oxy_numerator)$intensity ) /
                   (filter(tmp_squash, type != oxy_numerator)$intensity ))
    #}
  })
 # parallel::stopCluster(cl)
  data.frame(df_out)
}

volcano_plotter <- function(welch, phenotype_a, phenotype_b, title){
  out <- cbind(phenotype_a$mz,log2(apply(cbind(phenotype_a$fold,
                phenotype_b$fold), 1, max)),
                -log10(welch$p))
  out <- as.data.frame(out)
  colnames(out) <- c("mz","max_fold","p")
  out$diff <- "NO"
  out$diff[out$max_fold > 0.6 & out$p > -log10(0.05)] <- "UP"
  out$diff[out$max_fold < -0.6 & out$p > -log10(0.05)] <- "DOWN"
  as.data.frame(out)

  ggplot(data = out, aes(x = max_fold,
                                    y = p,
                                    col = diff)) +
    geom_point() +
    scale_color_manual(values = c("blue", "black", "red")) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    ggtitle("title")
}
