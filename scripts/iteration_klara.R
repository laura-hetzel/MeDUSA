#' @title Binning iterations
#' @description the code iterates the binprak function
#' over the data as long as the number of bins changes,
#' when the number of bins doesn't change anymore
#' iterations stop, maximal number of iterations is
#' set to 8
#' plot of the bin decrease can be created (optional)
#' @importFrom plyr join_all
#' @importFrom tidyr pivot_wider
#' @param df_list list of dataframes obtained from `binPeaks`
#' @param max_align (optional) value obtained from user input,
#' default value set to 8
#' @param bin_plot (optional) logical value deciding if plot
#' is created obtained from user input, default set to FALSE,
#' line plot shows reduction of bins per iteration
#' @export
iteration <- function(df_list, max_align = 8, bin_plot = F) {
  count <- 0L
  bins <- c()
  df <- plyr::join_all(df_list,
                       by = NULL,
                       type = "full", match = "all"
  )
  df <- df[, c("name", "mz", "intensity")]
  df <- pivot_wider(df,
                    names_from = "name",
                    id_cols = "mz",
                    values_from = "intensity"
  )

  df <- df[order(df$mz), ]
  bins <- c(bins, nrow(df))

  while (TRUE) {
    new <- binPeaks(df_list)
    count <- count + 1L
    new_df <- plyr::join_all(new,
                             by = NULL,
                             type = "full", match = "all"
    )
    new_df <- pivot_wider(new_df,
                          names_from = "name",
                          id_cols = "mz",
                          values_from = "intensity"
    )
    new_df <- new_df[order(new_df$mz), ]
    bins <- c(bins, nrow(new_df))
    if (nrow(df) == nrow(new_df)) {
      break
    }
    if (count == max_align) {
      break
    }
    df_list <- new
    df <- new_df
  }
  if (bin_plot == T) {
    df_bins <- as.data.frame(cbind(bins, 0:(count)))
    reduction_plot <- ggplot(df_bins) +
      geom_line(aes(x = V2, y = bins)) +
      labs(
        x = "Number of Alignments",
        y = "Number of bins"
      )
    df <- list(
      df = df,
      reduction_plot = reduction_plot
    )
    return(df)
  }
  return(df)
}
