#' @title Mass defect (MD) calculation
#' @description here a simplified MD, the mass' decimal number
#' @param dataframe data frame containing experimental data
#' @export
mass_defect_calculation <- function(dataframe) {
  MD <- dataframe$mz %% 1
  dataframe$MD <- MD
  return(dataframe)
}

#' @title Defining MD cut-off
#' @description the theoretical maximum MD of human metabolites
#' is used to calculate this cut off
#' @param filter_df obtained from the HMDB using `mass_defect_calculation`
#' @export
make_filter_list <- function(filter_df) {
  filter_df$mz <- 0.00112 * filter_df$mz + 0.01953
  return(filter_df)
}


#' @title Filtering the data by MD cut-off and inclusion list
#' @param MD_df obtained from from the HMDB using `mass_defect_calculation`
#' @param filter_df obtained from the `make_filter_list`
#' @param incl_list (optional) logical value from user input, default set
#' to FALSE, decides if inclusion list is used or not
#' @param incl_list_path (optional)path to inclusion list obtained from user input
#' or use of default inclusion list (McMillan et al., 2016)
#' @param mass_accuracy (optional) obtained from user input or use of
#' default value set to 0.01
#' @importFrom utils read.delim
MD_filter <- function(MD_df, filter_df, incl_list = FALSE,
                      incl_list_path = system.file("extdata/hmdb_inclusions_list_pos_McMillan.txt",
                      package = "sumR"), mass_accuracy = 0.01) {
  if (incl_list == T) {
    incl_list <- read.delim(file = incl_list_path, sep = "\t")
    incl_list_filtered <- incl_list[which(incl_list$inclusion == "y"), ]
    # compare compounds from inclusion list to the data
    cmp <- function(MD_df, incl_list_filtered, cutoff = mass_accuracy) {
      abs(incl_list_filtered - MD_df) <= cutoff
    }
    # assign comparable values a logical TRUE value
    match <- which(outer(incl_list_filtered$mz, MD_df$mz, cmp), arr.ind = TRUE)
    # return row names of TRUE rows
    filtered_rows <- rownames(MD_df)[as.numeric(levels(factor(match[, 2])))]
    # keep rows in MD_df if m/z defect is less than or equal to fitted value
    # OR value is in inclusion list
    filter_df <- MD_df[which((MD_df$MD <= filter_df$mz) | rownames(MD_df) %in% filtered_rows), ]
  }
  else {
    filter_df <- MD_df[which(MD_df$MD <= filter_df$mz)]
  }
  return(filter_df)
}

#' @title Plotting of the filtered data - m/z vs. MD
#' @param MD_df data frame obtained from the experimental data using `mass_defect_calculation`
#' @param MD_df_filtered data frame obtained from the experimental data using `MD_filter`
#' @export
plot_mz_MD <- function(MD_df, MD_df_filtered) {
  mz_removed <- nrow(MD_df) - nrow(MD_df_filtered)
  plot(MD_df_filtered$mz, MD_df_filtered$MD,
    cex.axis = 0.8,
    col = alpha("black", 0.5), pch = 20, cex = 0.8,
    ylim = c(0, 1), xlim = c(50, 1200), ylab = "MD", xlab = "m/z",
    main = "Filtered Data", sub = paste("datapoints removed = ", mz_removed),
    cex.lab = 0.8, cex.main = 0.8, cex.sub = 0.8
  )
}

#' @title Mass defect filter pipeline
#' @param dataframe containing the data, intensities and mz
#' @param mz_MD_plot logical variable deciding if plot is wanted per default set to TRUE
#' @importFrom stats na.omit
#' @importFrom utils read.delim
#' @export
MassDefectFilter <- function(dataframe, mz_MD_plot = TRUE) {
  # calculate the MD for all compounds
  md_df_exp <- mass_defect_calculation(dataframe)
  # make filtered list experimental data
  filter_list_exp <- make_filter_list(md_df_exp)
  # filter experimental data
  filtered_df_exp <- MD_filter(md_df_exp, filter_list_exp)
  # plotting
  if (mz_MD_plot == T) {
    plot_mz_MD(md_df_exp, filtered_df_exp)
  }}

