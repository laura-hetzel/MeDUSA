
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
#' @param metadata \cr
#'   DataFrame : Metadata-Obj of Samples
#' @param annotation \cr
#'   String: Colname of metadata to annotate on.
#' @param title
#'   String: Plot title
#' @param save_file
#'   String: Path to save the plot.
#'   TRUE: Save to default output directory
#'   NA: Do not save
#' @param cluster
#'   c("row"=Boolean,"col"=Boolean): Cluster data? Defaults to False
#'
#' @export
plot_heatmap <- function(input_mz_obj, metadata = NULL, annotation = "phenotype",
                          title = "HeatMap of intensities", save_file = TRUE, cluster = NULL, ...){

  plot_obj <- dplyr::select(input_mz_obj, -mz)
  row.names(plot_obj) <- input_mz_obj$mz

  if (!is.null(annotation) && !is.null(metadata)){
    if (sum( sort(metadata$sample_name) !=
             sort(colnames(plot_obj)))) {
      stop("ERROR: plot_heatmap: annotation sample_name does not match mz_obj data")
    }
    annotation_int = data.frame(row.names = metadata$sample_name, metadata[annotation])
  } else {
    annotation_int <- NULL
  }
  if (!is.na(save_file) && save_file == T){
    save_file = paste(local.output_dir,local.dir_sep(),title,".png",sep = "")
  }
  if (is.null(cluster)){
    cluster <- c("row"=F,"col"=F)
  }

  pheatmap::pheatmap( t(plot_obj),
                      cluster_rows = cluster["row"],
                      cluster_cols = cluster["col"],
                      annotation_row = annotation_int,
                      main = title,
                      filename = save_file,
                      ...)
}
