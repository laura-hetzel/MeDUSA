#' @title Mass defect (MD) calculation
#' @description here a simplified MD, decimal number of the mass
#' @param dataframe data frame containing experimental data
#' df must contain an mz column
mass_defect_calculation <- function(dataframe) {
  MD <- dataframe$mz %% 1
  dataframe$MD <- MD
  return(dataframe)
}

#' @title Defining MD cut-off
#' @description the theoretical maximum MD of human metabolites
#' is used to calculate this cut off
#' @param filter_df obtained from the HMDB using `mass_defect_calculation`
make_filter_list <- function(filter_df) {
  filter_df$mz <- 0.00112 * filter_df$mz + 0.01953
  return(filter_df)
}


#' @title Filtering the data by MD cut-off,(optional) and inclusion list
#' @param MD_df obtained from from the HMDB using `mass_defect_calculation`
#' @param filter_df obtained from the `make_filter_list`
#' @param incl_list (optional) logical value from user input, default set
#' to FALSE, decides if inclusion list is used or not
#' @param incl_list_path (optional, use only if incl_list = T)path to
#' inclusion list obtained from user input, default inclusion list from
#' @param mass_accuracy (optional, use only if incl_list = T) obtained
#' from user input, default value set to 0.01
#' @importFrom utils read.delim
MD_filter <- function(MD_df, filter_df, incl_list,
                      incl_list_path = system.file("extdata/hmdb_inclusions_list_pos_McMillan.txt",
                                                                                package = "sumR"),
                      mass_accuracy) {
  if (incl_list == T) {
    list <- read.delim(file = system.file("extdata/hmdb_inclusions_list_pos_McMillan.txt",
                                          package = "sumR"), sep = "\t")
    incl_list_filtered <- list[which(list$inclusion == "y"), ]
    # compare compounds from inclusion list to the data
    cmp <- function(MD_df, incl_list_filtered, cutoff = mass_accuracy) {
      abs(incl_list_filtered - MD_df) <= cutoff
    }
    # assign comparable values a logical TRUE value
    match <- which(outer(incl_list_filtered$mz, MD_df$mz, cmp), arr.ind = TRUE)
    # return row names of TRUE rows
    filtered_rows <- rownames(MD_df)[as.numeric(levels(factor(match[, 2])))]
    ## keep rows in MD_df if m/z defect is less than or equal to fitted value
    ## OR value is in inclusion list
    filtered_df <- MD_df[which((MD_df$MD <= filter_df$mz) | rownames(MD_df) %in% filtered_rows), ]
  }
  else {
    filtered_df <- MD_df[which(MD_df$MD <= filter_df$mz), ]
  }
  return(filtered_df)
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
#' @description uses the mass defect functions to filter
#' the data with one command
#' @param dataframe containing the data, intensities and mz
#' @param mz_MD_plot logical variable deciding if plot is
#' created per default set to TRUE
#' @importFrom stats na.omit
#' @importFrom utils read.delim
#' @export
massDefectFilter <- function(exp, plot = FALSE, incl_list = FALSE,
                             incl_list_path = system.file("extdata/hmdb_inclusions_list_pos_McMillan.txt",
                                                          package = "sumR"),
                             mass_accuracy = 0.01 ) {
  if (!validateExperiment(exp)) return(NULL)
  dataframe <- as.data.frame(rowData(exp))
  # calculate the MD for all compounds
  md_df <- mass_defect_calculation(dataframe)
  # make filtered list experimental data
  filter_list <- make_filter_list(md_df)
  # filter experimental data
  filtered_df <- MD_filter(MD_df = md_df, filter_df = filter_list,
                               incl_list = incl_list, mass_accuracy = mass_accuracy,
                               incl_list_path = incl_list_path)
  # plotting
  if (plot) {
    plot_mz_MD(md_df, filtered_df)
  }

  exp <- exp[rownames(filtered_df), ]
  rowData(exp) <- DataFrame(filtered_df)
  exp
}

