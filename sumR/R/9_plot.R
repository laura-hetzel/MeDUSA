
# *** Volcano Plot -----------------------------------------------------
#' Welch & FoldChange Volcano Plot
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
#' @export
plot_volcano <- function(welch, fold_change, title = "Volcano Plot"){
  if ( welch$mz != fold$mz){
    stop("ERROR: plot_volcano : Welch mz does not match FoldChange mz")
  }
  df <- cbind(welch$mz, log2(welch$p), log2(fold_change$fold))

  colnames(df) <- c("mz","p","fold")
  df$diff <- "NO"
  df$diff[out$max_fold > 0.6 & df$p > -log10(0.05)] <- "UP"
  df$diff[out$max_fold < -0.6 & df$p > -log10(0.05)] <- "DOWN"
  as.data.frame(df)

  ggplot( data = df,
          aes(x = fold, y = p, col = diff)) +
        geom_point() +
          scale_color_manual(values = c("blue", "black", "red")) +
          geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
          geom_hline(yintercept = -log10(0.05), col = "red") +
          ggtitle(title)

  local.save_plot(paste("Volcano",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
}

# *** Volcano Magic -----------------------------------------------------
#' MZLOG-OBJ Volcano Magic
#' Create Welch, Fold & Volcano plot, with default inputs
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : Log2 of Input MZ-Obj
#' @param phenotype_a \cr
#'   DataFrame: Filtered metadata for phenotypeA
#' @param phenotype_b \cr
#'   DataFrame: Filtered metadata for phenotypeB
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
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
#' @export
mz_analysis_volcano_magic <- function(input_mzlog_obj, phenotype_a, phenotype_b, cores = 2 ){
  welch <- mzlog_analysis_welch(input_mzlog_obj, phenotype_a, phenotype_b, cores)
  fold  <- mzlog_analysis_fold(input_mzlog_obj, phenotype_a, phenotype_b, cores)

  plot_volcano(welch, fold)

}

# *** Heat Map  -----------------------------------------------------
#' MZLOG-OBJ Plot Heat map
#'
#' Create Welch, Fold & Volcano plot, with default inputs
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param annotation \cr
#'   DataFrame: Two columns from metadata: (sample_name ; [annotation_column]) * optional
#' @param title
#'   String: Plot title
#' @param save_file
#'   String: Path to save the plot. Defaults to default output directory
#' @param cluster
#'   c("row"=Boolean,"col"=Boolean): Cluster data? Defaults to False
#'
#' @examples
#' Create a HeatMap noting Phenotype
#'   plot_heatmap(mz_object, select(metadata, "sample_name", "phenotype"))
#'
#'
#' @export
plot_heatmap <- function(input_mz_obj, annotation = NULL, title = "HeatMap of intensities", save_file = NULL, cluster = NULL...){
  if (!is.null(annotation)){
    if (sort(annotation$sample_name) != sort(dplyr::select(input_mz_obj,-mz))){
      stop("ERROR: plot_heatmap: annotation sample_name does not match mz_obj data")
    }
    if (dim(annotation)[2] != 2){
      stop(paste("ERROR: plot_heatmap: annotation parameter expected 2 columns. Got: ", dim(annoation[2]),dir_sep=""))
    }
    annotation = data.frame(row.names = annotation$sample_name, dplyr::select(annotation, -sample_name)))
  }
  if (is.null(save_file)){
    save_file = paste(local.output_dir,local.dir_sep(),title,".png",sep = "")
  }
  if (is.null(cluster)){
    cluster <- c("row"=F,"col"=F)
  }

  pheatmap::pheatmap( t(input_mz_obj),
                      cluster_rows = cluster["row"],
                      cluster_cols = cluster["col"],
                      annotation_row = annotation,
                      main = title,
                      filename = save_file,
                      ...)
}
