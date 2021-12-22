#' @title Retrieve HMDB data in a dataframe format
#' @param file (optional) String location of the HMDB file
#' @importFrom stats na.omit
#' @importFrom tidyverse read_delim
#' @export
get_hmdb_file <- function(file = NULL){
  if (is.null(file)) {
    # No file given, must download from HMDB
    # TODO: Find out how to download from HMDB
    url <- "https://hmdb.ca/metabolites.csv?action=index&c=hmdb_id&controller=metabolites&d=up&endogenous=1&filter=true&utf8=%E2%9C%93"
  } else {
    #' Downloading the data and converting it to a data frame
    dataframe <- read_delim(file = hmdb_data, delim = ",")            #82.595 rows as should be
    #removing rows with no Mass data
    dataframe <- na.omit(dataframe)
  }
  dataframe
}


#' @title Mass defect filter function
#' @param dataframe Dataframe obtained from HMDB using `get_hmdb_file`
#' @importFrom dplyr mutate
#' @export
mdf_function <- function(dataframe){
  mz_col <- dataframe$MONO_MASS
  MD <- mz_col %% 1
  #to add the MD to the new data frame
  MD_df <- dplyr::mutate(dataframe, MD)
  MD_df <- MD_df[,7:8]
  MD_df
}




#' @title
#' @param dataframe MD_df obtained from the HMDB using 'mdf_function' 
#' @export
create_fit_df <- function(MD_df){
  fit <- as.data.frame(apply(MD_df[,1], 1, function(row){
    y <- 0.00112 * row + 0.01953
    return(y)
  }))

  colnames(fit) <- "mz"
  fit
}

#' @title
#' @param dataframe MD_df obtained from the HMDB using 'mdf_function' or filtered_MD_df obtained from the HMDB using 'filter_function'
#' @export
plot_m/z_MD <- function(MD_df, title){
  plot(MD_df$MONO_MASS, MD_df$MD, cex.axis = 0.8,
       col = alpha("black", 0.5), pch = 20, cex = 0.8,
       ylim = c(0,1), xlim = c(50,1200), ylab = "MD", xlab = "m/z",
       main = title, cex.lab = 0.8, cex.main = 0.8)
}


inclusion_list_function <- function(file = NULL){
  # Inclusion list creation
  # default file from McMillan paper 
  # no file given must use this one 
  if (is.null(file)) {
    in_list_data <- read.table("hmdb_inclusions_list_pos_McMillan.txt",  header = T, check.names = F,
                               row.names = 1, sep = "\t")
  } else {
    in_list_data <- read.table(file,  header = T, check.names = F,
                               row.names = 1, sep = "\t")
  }
  #' set mass accuracy (in Daltons) for inclusion list matching
  ma <- 0.01



  #find masses in inclusion list in MD_df
  in_list <- hmdb_list[which(hmdb_list$"inclusion" == "y"),]                  #if in the input data inclusion list is put to y (yes) it is included here
  cmp <- function(MD_df, in_list, cutoff=ma){abs(in_list - MD_df) <= cutoff}  #compare compounds from from inclusion list
  match <- which(outer(in_list$mz, MD_df$MONO_MASS, cmp), arr.ind = TRUE)     # assign comparable values a logical TRUE value
  fc <- factor(match[,2])                                                     # make a factor out of the column 2
  I <- rownames(MD_df)[as.numeric(levels(fc))]                               # return row names of TRUE rows

  #keep rows in MD_df if m/z defect is less than or equal to fitted value OR value is in inclusion list

  filtered <- MD_df[which((MD_df$"MD" <= fit$"mz") | rownames(MD_df) %in% I),]
  filtered
}


pipeline <- function(file = NULL){
  dataframe <- get_hmdb_file(file)
  md_df <- mdf_function(dataframe)
  fit <- create_fit_df(md_df)
  plot_mono_mass(md_df)
  filtered_df <- filter_function(md_df)
  plot_mono_mass(filtered_df)
}
