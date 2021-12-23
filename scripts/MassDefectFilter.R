#' @title Retrieve HMDB data in a data frame format - default = endogenous metabolites 
#' @param file (optional) String location of the HMDB file
#' @importFrom stats na.omit
#' @importFrom tidyverse read_delim
#' @export
get_hmdb_file <- function(file){
  #if (is.null(file)) {
    # No file given, must download from HMDB
    # TODO: Find out how to download from HMDB
    #url <- "https://hmdb.ca/metabolites.csv?action=index&c=hmdb_id&controller=metabolites&d=up&endogenous=1&filter=true&utf8=%E2%9C%93"
  #} else {
    #' Downloading the data and converting it to a data frame
    dataframe <- read_delim("end.metabolites.txt", delim = ",")            #82.595 rows as should be
    #removing rows with no Mass data
    dataframe <- na.omit(dataframe)
  #}
  return(dataframe)
}

#' @title Import experimental data
#' @param file String location of the experimetal file
#' @importFrom stats na.omit
#' @importFrom tidyverse read_delim
#' @export
get_data <- function(file){
  #' Downloading the experimental data and converting it to a data frame
  ex_data_df <- read_delim("889C0121170_NS_005_NEG_7_cvd_d5_pm.txt", col_names = F)
  colnames(ex_data_df) <- c("p", "MONO_MASS", "l", "Intesity","e")
  return(ex_data_df)
}

#' @title Mass defect calculation function
#' @param dataframe Data frame obtained from HMDB using `get_hmdb_file`
#' @importFrom dplyr mutate
#' @export
mass_defect_calculation <- function(dataframe){
  mz_col <- dataframe$MONO_MASS
  MD <- mz_col %% 1
  #to add the MD to the new data frame
  MD_df <- dplyr::mutate(dataframe, MD)
  return(MD_df)
}

#' @title Defining Mass Defect filter 
#' @description  the theoretical maximum Mass Defect of human metabolites is used to calculate this cut off
#' @param dataframe MD_df obtained from the HMDB using `mass_defect_function`
#' @export
make_filter_list <- function(MD_df){
  filter_list <- as.data.frame(apply(MD_df[,"MONO_MASS"], 1, function(row){
    y <- 0.00112 * row + 0.01953
    return(y)
  }))
  colnames(filter_list) <- "mz"
  return(filter_list)
}

#' @title Retrieve inclusion list in a data frame format (default or user input)
#' @param file (optional) String location of the inclusion list
#' @importFrom tidyverse read_delim
#' @export
get_inclusion_list <- function(file){
  # Inclusion list creation
  # default file from McMillan paper 
  # no file given must use this one 
 # if (is.null(file)) {
    incl_list <- read_delim("hmdb_inclusions_list_pos_McMillan.txt", delim = "\t")
  #} else {
   # incl_list <- read_delim(file, delim = "\t")
  #}
  return(incl_list)
}

#' @title Filtering the data by  mass defect 
#' @param dataframe MD_df obtained from from the HMDB using `mass_defect_function`
#' @param dataframe incl_list obtained from the `get_inclusion_list`
#' @param dataframe filter_list obtained from the `create_fit_df`
#' @importFrom dplyr filter
MD_filter <- function(MD_df, filter_list, incl_list){
  #' set mass accuracy (in Daltons) for inclusion list matching
  ma <- 0.01
  #find masses in inclusion list in MD_df
  incl_list_filtered <- incl_list[which(incl_list$"inclusion" == "y"),]                  #if in the input data inclusion list is put to yes it is included here
  cmp <- function(MD_df, incl_list_filtered, cutoff = ma){abs(incl_list_filtered - MD_df) <= cutoff}  #compare compounds from inclusion list to the data
  match <- which(outer(incl_list_filtered$mz, MD_df$MONO_MASS, cmp), arr.ind = TRUE)     # assign comparable values a logical TRUE value
  filtered_rows <- rownames(MD_df)[as.numeric(levels(factor(match[,2])))]                               # return row names of TRUE rows
  #keep rows in MD_df if m/z defect is less than or equal to fitted value OR value is in inclusion list
  filtered_df <- MD_df[which((MD_df$MD <= filter_list$mz) | rownames(MD_df) %in% filtered_rows),]
  return(filtered_df)
}

#' @title Plotting of the data - m/z vs. MD 
#' @param dataframe MD_df obtained from the HMDB using `mass_defect_function` or
#' @param dataframe filtered_MD_df obtained from the HMDB using `MD_filter`
#' @export
plot_mz_MD <- function(MD_df, title){
  plot(MD_df$MONO_MASS, MD_df$MD, cex.axis = 0.8,
       col = alpha("black", 0.5), pch = 20, cex = 0.8,
       ylim = c(0,1), xlim = c(50,1200), ylab = "MD", xlab = "m/z",
       main = title, cex.lab = 0.8, cex.main = 0.8)
}

pipeline <- function(file = NULL){
  #import hmdb
  dataframe_hmdb <- get_hmdb_file(file)
  #calculate MD HMDB
  md_df_hmdb <- mass_defect_calculation(dataframe_hmdb)
  #make filtered list hmdb
  filtered_list_hmdb <- make_filter_list(md_df_hmdb)
  #import experimental data
  ex_data_df <- get_data(file)
  # calculate MD experimental data 
  md_df_exp <- mass_defect_calculation(ex_data_df)
  #make filtered list experimental data
  filtered_list_exp <- make_filter_list(md_df_exp)
  # import inclusion list
  in_list <- get_inclusion_list(file)
  # filter hmdb
  filtered_df_hmdb <- MD_filter(md_df_hmdb, filtered_list_hmdb, in_list)
  # filter experimental data 
  filtered_df_exp <- MD_filter(md_df_exp, filtered_list_exp, in_list)
  # plotting
  plot_mz_MD(md_df_exp, title = "Raw")
  plot_mz_MD(filtered_df_exp, title = "Filtered")
}

pipeline()
