#' @title HMDB - loading and plotting
#' @description This function loads the hmdb xml file directly over HTTP, plots MZ vs. MD and creats the cutoff for the MD filter.
#' @usage
#' @importFrom tidyverse read_delim
#' @export
#' @example
#' @author Amar Poeran & Klara Blaschitz

#' Download required packages
library(tidyverse)

#' Input by user
#' download the data to your wd and enter file name here
hmdb_data <- "end.metabolites.txt"       #downloaded txt file with only endogenous metabolites

#' inclusion list read in
in_list_data <- "hmdb_inclusions_list_pos_McMillan.txt"

#' set mass accuracy (in Daltons) for inclusion list matching
ma <- 0.01


#' Downloading the data and converting it to a data frame
dataframe <- read_delim(file = hmdb_data, delim = ",")            #82.595 rows as should be
#removing rows with no Mass data
dataframe <- na.omit(dataframe)                         #82.580, 15 rows removed because Mass was NA

#' Function to calculate Mass Defect for filtering
mz_col <- dataframe$MONO_MASS
MD <- function(dataframe, mz_col){
  MD <- mz_col %% 1
  #to add the MD to the new data frame
  MD_df <- dplyr::mutate(dataframe, MD)
}
#' Apply the function
MD_df <- MD(dataframe, mz_col)
MD_df <- MD_df[,7:8]

#' integrating the cutoff line in the data
f <- function(x){
  y <- 0.00112 * x + 0.01953    # cut-off linear equation
  return(y)
}

#MD_m <- matrix(MD_df[,7])

fit <- as.data.frame(apply(MD_df[,1],1,f))
colnames(fit) <- "mz"

#' Inclusion list creation
#' reproduced from the MacMillan salt cluster removal paper

hmdb_list <- read.table(in_list_data,  header = T, check.names = F,
              row.names = 1, sep = "\t")


#find masses in inclusion list in MD_df
in_list <- hmdb_list[which(hmdb_list$"inclusion" == "y"),]                  #if in the input data inclusion list is put to y (yes) it is included here
cmp <- function(MD_df, in_list, cutoff=ma){abs(in_list - MD_df) <= cutoff}  #compare compounds from from inclusion list
match <- which(outer(in_list$mz, MD_df$MONO_MASS, cmp), arr.ind = TRUE)     # assign comparable values a logical TRUE value
fc <- factor(match[,2])                                                     # make a factor out of the column 2
I <- rownames(MD_df)[as.numeric(levels(fc))]                               # return row names of TRUE rows

#keep rows in MD_df if m/z defect is less than or equal to fitted value OR value is in inclusion list

filtered <- MD_df[which((MD_df$"MD" <= fit$"mz") | rownames(MD_df) %in% I),]

#' plotting results
plot(MD_df$"MONO_MASS", MD_df$"MD", cex.axis = 0.8,
     col = alpha("black", 0.5), pch = 20, cex = 0.8,
     ylim = c(0,1), xlim = c(50,1200), ylab = "MD", xlab = "m/z",
     main = "Raw", cex.lab = 0.8, cex.main = 0.8)


plot(filtered$"MONO_MASS", filtered$"MD", cex.axis = 0.8,
     col = alpha("black", 0.5), pch = 20, cex = 0.8,
     ylim = c(0,1),xlim = c(50,1200), ylab = "MD", xlab = "m/z",
     main = "Filtered", cex.lab = 0.8, cex.main = 0.8)


