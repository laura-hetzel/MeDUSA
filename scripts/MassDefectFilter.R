#' @title Retrieve HMDB data in a data frame format - default = endogenous metabolites 
#' @param file (optional) String location of the HMDB file
#' @importFrom stats na.omit
#' @importFrom readr read_delim
#' @export
get_hmdb_file <- function(file){
  #if (is.null(file)) {
    # No file given, must download from HMDB
    # TODO: Find out how to download from HMDB
    #url <- "https://hmdb.ca/metabolites.csv?action=index&c=hmdb_id&controller=metabolites&d=up&endogenous=1&filter=true&utf8=%E2%9C%93"
  #} else {
    #' Downloading the data and converting it to a data frame
    #dataframe <- read.delim("end.metabolites.txt", sep = ",")            
    dataframe <- read_delim("end.metabolites.txt", delim = ",")
    colnames(dataframe) <- c("HMDB_ID", "Name", "SMILES", "INCHIKEY", "Chemical_Formula", "Average_Mass", "mz")
    #removing rows with no Mass data
    dataframe <- na.omit(dataframe)
  #}
  return(dataframe)
}

#' @title Import experimental data
#' @param file String location of the experimental file
#' @importFrom readr read_delim
#' @export
get_data <- function(file){
  #' Downloading the experimental data and converting it to a data frame
  #ex_data_df <- read.delim("peaks.csv", sep = ",")
  ex_data_df <- read_delim("peaks.csv", col_names = T)
  return(ex_data_df)
}

#' @title Mass defect calculation function
#' @param dataframe Data frame obtained from HMDB using `get_hmdb_file`
#' @export
mass_defect_calculation <- function(dataframe){
  MD <- dataframe$mz %% 1
  #to add the MD to the new data frame
  dataframe$MD <- MD
  return(dataframe)
}



#' @title Defining Mass Defect filter 
#' @description  the theoretical maximum Mass Defect of human metabolites is used to calculate this cut off
#' @param filtered_df obtained from the HMDB using `mass_defect_calculation`
#' @export
make_filter_list <- function(filtered_df){
  filtered_df$mz <- 0.00112 * filtered_df$mz + 0.01953
  return(filtered_df)
}

#' @title Filtering the data by  mass defect 
#' @param MD_df obtained from from the HMDB using `mass_defect_calculation`
#' @param filtered_df obtained from the `make_filter_list`
#' @param incl_list path to inclusion list obtained from user input or use of default inclusion list (McMillan et al., 2016)
#' @param mass_accuracy obtained from user input or use of default value set to 0.01
#' @importFrom utils read.delim
MD_filter <- function(MD_df, filtered_df, 
                      incl_list = "hmdb_inclusions_list_pos_McMillan.txt", mass_accuracy = 0.01){
  incl_list <- read.delim(file = incl_list, sep = "\t")
  #find masses in inclusion list in MD_df
  incl_list_filtered <- incl_list[which(incl_list$inclusion == "y"),]                  #if in the input data inclusion list is put to yes it is included here
  cmp <- function(MD_df, incl_list_filtered, cutoff = mass_accuracy){abs(incl_list_filtered - MD_df) <= cutoff}  #compare compounds from inclusion list to the data
  match <- which(outer(incl_list_filtered$mz, MD_df$mz, cmp), arr.ind = TRUE)     # assign comparable values a logical TRUE value
  filtered_rows <- rownames(MD_df)[as.numeric(levels(factor(match[,2])))]                               # return row names of TRUE rows
  #keep rows in MD_df if m/z defect is less than or equal to fitted value OR value is in inclusion list
  filtered_df <- MD_df[which((MD_df$MD <= filtered_df$mz) | rownames(MD_df) %in% filtered_rows),]
  return(filtered_df)
}

#' @title Plotting of the data - m/z vs. MD 
#' @param MD_df dataframe obtained from the experimental data using `mass_defect_calculation`
#' @param filtered_MD_df dataframe obtained from the experimental data using `MD_filter`
#' @export 
plot_mz_MD <- function(MD_df_filtered, MD_df){
  mz_removed <- nrow(MD_df) - nrow(MD_df_filtered)
  plot(MD_df_filtered$mz, MD_df_filtered$MD, cex.axis = 0.8,
       col = alpha("black", 0.5), pch = 20, cex = 0.8,
       ylim = c(0,1), xlim = c(50,1200), ylab = "MD", xlab = "m/z",
       main = "Filtered Data", sub = paste("datapoints removed = ", mz_removed),
       cex.lab = 0.8, cex.main = 0.8, cex.sub = 0.8)
}

#' @title Export results in a csv file 
#' @param dataframe filtered_MD_df obtained from `MD_filter`
result_output <- function(filtered_df){
  #write filtered_df in csv 
  write.csv(filtered_df,"filtered_data.csv", row.names = F)
}


#' @title Mass defect filter pipeline 
#' @param file path to file that contains the input data 
#' @param mz_MD_plot logical variable deciding if plot is wanted per default set to TRUE 
#' @importFrom stats na.omit
#' @importFrom readr read_delim
#' @importFrom utils read.delim
#' @export
MassDefectFilter <- function(file = NULL, mz_MD_plot = TRUE){
  ex_data_df <- get_data(file)
  # calculate the MD for all compounds
  md_df_exp <- mass_defect_calculation(ex_data_df)
  # make filtered list experimental data
  filtered_list_exp <- make_filter_list(md_df_exp)
  # filter experimental data 
  filtered_df_exp <- MD_filter(md_df_exp, filtered_list_exp)
  # plotting
  if (mz_MD_plot == T) {
    plot_mz_MD(filtered_df_exp, md_df_exp)
  }
}

MassDefectFilter()





