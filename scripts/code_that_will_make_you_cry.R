ppm_calc <- function (mass1, mass2) {
  #' ppm_calc calculated the parts per million error between two different masses
  #' 
  #' @examples 
  #' ppm_calc(mass1, mass2)
  #' 
  #' @export
  ppm_error <- ((mass1 - mass2)/mass1) * 1e6
  return(ppm_error)
}
align_check <- function(data_frame_fn, xcoords = c(-20, 0)) {
  #' align_check makes sure that all the m/z values are aligned/binned correctly
  #' 
  #' align_check takes 2 arguments (1 optional), a dataframe of peaks with mz column,
  #' and an optinal argument for the coordinates of the plot to be zoomed in on
  #' 
  #' align_check outputs a list of three elements:
  #' 1- Dataframe of 1 column containing the ppm error values
  #' 2- boxplot of the ppm erro values with xcoords zoomed in to -20,0 (default)
  #' 3- table of summary stats of ppm error values
  #' 
  #' @importFrom ppm_calc
  #' @examples 
  #' 
  #' 
  #' @export
  odd_ind_fn <- seq(3, length(data_frame_fn$mz), 2)
  even_ind_fn <- seq(2, length(data_frame_fn$mz),2)
  ppm_err_fn <- ppm_calc(data_frame_fn$mz[even_ind_fn], data_frame_fn$mz[odd_ind_fn]) %>% 
    as_tibble_col( column_name = "ppm_error") 
  ppm_err_plot_fn <- ggplot(ppm_err_fn, aes(x = ppm_error)) +
    geom_boxplot() +
    coord_cartesian(xlim = xcoords) +
    ggtitle("ppm error boxplot") +
    theme_classic(base_size = 20)
  ppm_err_summary_fn <- summary(ppm_err_fn)
  ppm_err_list <-  list(ppm_error_df = ppm_err_fn,
                        boxplot = ppm_err_plot_fn, 
                        summary_stats =  ppm_err_summary_fn
  )
  #paste(nrow(data_frame) - ) add a print
  return(ppm_err_list)
}

#'BG_filter is used to conduct background removal (done by amar)
#'
#'BG_filter is a function used to filter a dataframe of m/z values and intensities. 
#'This function allows the user to remove m/z peaks of measured samples, if they correspond to the same m/z of a blank within
#'a specified range (for example within 150% of the blank intensity).
#'
#'@param dataframe
#'@param mz_column
#'
#'
#'@importFrom
#'
#'
#'@examples
#'
#'
#'@export

BG_filter<- function(dataframe, limit=2.5, blank_regex ="blank", mz_regex="mz", sample_regex="^d|^ec|qc"){
  blank_cols <- grep(blank_regex, names(dataframe),ignore.case = TRUE)
  
  sample_cols <- grep(sample_regex, names(dataframe),ignore.case = TRUE)
  dataframe$background <- rowMeans(dataframe[, blank_cols], na.rm = TRUE)
  #Compare sample_cols with the mean, replace them by 0 if they are below thresh
  #From what I understand up until now: Sweep creates the Intensity Ratio as a summary statistic, which sweeps across the sample columns
  # The function acrosse those columns is basically a division by the corresponding Background value.
  #I'm gettin errors related to: Dimensions of matrix need to be equal to the function input data
  #https://stackoverflow.com/questions/3444889/how-to-use-the-sweep-function
  dataframe[sample_cols][sweep(dataframe[sample_cols], 1, dataframe$background, `/`) <= limit] <- 0
  result <- dataframe[, -blank_cols]
  return(result)
  #ahmeds notes: Need to review this later, maybe replace swap with map_dfc function
}



#' MD_filter for data filtering based on Mass Defect
#' 
#'
#' 
#' @importFrom dplyr mutate select
#' 
#' @param dataframe is the dataframe on which KM calculations need to be conducted
#' @param mz_col is the mass to charge column within the specified dataframe
#' @param a is the coefficient of the linear equation used for data filtering
#' @param b is the addition value of the linear equation used for data filtering, 
#' both a and b can be calculated using the linear_equation() function included in this package.
#' 
#' @examples 
#' MD_filter(master_df_pos, master_df_pos$mz)
#' 
#' @export

# Use Only this for now ---------------------------------------------------

MD_filter<- function(dataframe, mz_col, a = 0.00112, b = 0.01953){ #SOlves problem of customizable lin. eq as well: In case HMDB updates data
  #In-function MD calculation 
  MZ<- mz_col
  MZR<- trunc(mz_col, digits = 0)#Either floor() or trunc() can be used for this part.
  MD<- as.numeric(MZ-MZR)
  MD.limit<- b + a*mz_col
  dataframe<- dataframe%>%
    dplyr::mutate(MD, MZ, MD.limit)%>%
    dplyr::select(MD, MZ, MD.limit)
  #This approach uses the linear equation to create an additional column and based
  #on this column, we will be filtering out data.
  
  #Below are in total 3 diffferent plots: 1 of the starting dataframe, the second of the filtered dataframe
  #the third is basically added to the first, in order to highlight the "to be removed" datapoints as a means of tagging them prior to removal.
  
  highlight_df <- dataframe %>% filter(MD >= MD.limit) #Notice how this is the exact opposite from the 
  #"filtered" dataframe below. This one filters everything ABOVE the limit, as everything below will be kept and we want 
  #those datapoints which will be removed, highlighted. Not those which will remain.
  
  #DO the same for inlcusion list of HMDB, We need to make: 1, Linear eq. 2, inclusion list plot 3, metaboshiny. 4, R package.
  
  MD_plot<- ggplot(data=dataframe, aes(x=MZ, y=MD))+
    geom_point()+
    geom_point(data=highlight_df, aes(x=MZ,y=MD), color='red')+#I added this one, so the data which will be removed will be highlighted in red.
    ggtitle(paste("Unfiltered MD data"))# deparse(substitute(dataframe)) if you want to add the dataframe name, acc to St. Ovflw
  #stat_smooth(method="lm", se=FALSE)-> For linear line through the plot, but may not be necessary to show
  
  
  #Creating a filtered dataframe:
  filtered<- dataframe%>%
    filter(MD <= MD.limit)# As I understood: Basically all are coordinates. The maxima equation basically gives coordinates 
  #for the m/z values (x and MD = y). If it exceeds the equivalent coordinate of Y (which is MD) for the linear equation, it will be filtered.
  
  
  
  MD_plot_2<- ggplot(data=filtered, aes(x=MZ, y=MD))+ #Filtered is basically the second dataframe, #which subsets datapoints with an Y value (which is the MD), below the linear equation MD...
    geom_point()+
    ggtitle(paste("Filtered MD data"))
  #stat_smooth(method="lm", se=FALSE) -> For linear line through the plot, but may not be necessary to show
  N_Removed_datapoints <- nrow(dataframe) - nrow(filtered)#To determine the number of peaks removed
  print(paste("Number of peaks removed:", N_Removed_datapoints))
  MD_PLOTS<-list(preMD_df = dataframe, postMD_df = filtered, preMD_plot =  MD_plot, postMD_plot = MD_plot_2)
  return(MD_PLOTS)
}

FUN_text_exporter <- function(rnamed_cnamed_df,file_path,extension) {
  #takes three arguments,
  #df: with rownamed mz,colnamed samples
  #path : desire output path ,should be a string
  #custom_extension : desired custom extension to files (aside from .txt), string
  #takes teh df, separates it into mz,int files for each sample (1st loop)
  #export each file with desired file path,and custom extension (2nd loop)
  #UNALIGNED ! UNLIMITED POWAAAAAH
  all_spectra <- list()
  c_names <- colnames(rnamed_cnamed_df)
  r_names <- row.names(rnamed_cnamed_df)
  message("column /row names extracted , commencing 1st boop de loop")
  for (i in seq_along(c_names)) {
    assign(paste0(c_names[i],'_BG'), 
           as.data.frame(cbind(mz = as.numeric(r_names), 
                               int = rnamed_cnamed_df[ , c_names[i]])))
  }
  message("first loop successiful, commecning 2nd boop de loop")
  
  all_spectra <- lapply(ls()[grepl(pattern = '_BG', ls())], get, envir = environment())
  message(paste("Length of file list: ",length(all_spectra)))
  
  filenames <- c_names
  for (i in seq_along(filenames)) {
    write.table(subset(all_spectra[[i]], all_spectra[[i]][2] > 0),
                file = paste0(file_path,filenames[i],
                              extension,'.txt'),
                col.names = FALSE, row.names = FALSE, sep = '\t')
  }
  message("everything is done")
}




condition <- function (mass, intensities, samples, tolerance) 
{
  if (anyDuplicated(samples)) {
    return(NA)
  }
  meanMass <- mean(mass)
  if (any(abs(mass - meanMass)/meanMass > tolerance)) {
    return(NA)
  }
  meanMass
}
#Binning stuff ------------------------
##binning dependency 1

condition <- function (mass, intensities, samples, tolerance) 
{
  if (anyDuplicated(samples)) {
    return(NA)
  }
  meanMass <- mean(mass)
  if (any(abs(mass - meanMass)/meanMass > tolerance)) {
    return(NA)
  }
  meanMass
}

##Binning dependency 2

binning <- function (mass, intensities, samples, tolerance) 
{
  n <- length(mass)
  d <- diff(mass)
  nBoundaries <- max(20L, floor(3L * log(n)))
  boundary <- list(left = double(nBoundaries), right = double(nBoundaries))
  currentBoundary <- 1L
  boundary$left[currentBoundary] <- 1L
  boundary$right[currentBoundary] <- n
  while (currentBoundary > 0L) {
    left <- boundary$left[currentBoundary]
    right <- boundary$right[currentBoundary]
    currentBoundary <- currentBoundary - 1L
    gaps <- d[left:(right - 1L)]
    gapIdx <- which.max(gaps) + left - 1L
    l <- condition(mass = mass[left:gapIdx], intensities = intensities[left:gapIdx], samples = samples[left:gapIdx], tolerance = tolerance)
    if (is.na(l[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    }
    else {
      mass[left:gapIdx] <- l
    }
    r <- condition(mass = mass[(gapIdx + 1L):right], intensities = intensities[(gapIdx + 1L):right], samples = samples[(gapIdx + 1L):right], tolerance = tolerance)
    if (is.na(r[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- gapIdx + 1L
      boundary$right[currentBoundary] <- right
    }
    else {
      mass[(gapIdx + 1L):right] <- r
    }
    if (currentBoundary == nBoundaries) {
      nBoundaries <- floor(nBoundaries * 1.5)
      boundary$left <- c(boundary$left, double(nBoundaries -currentBoundary))
      boundary$right <- c(boundary$right, double(nBoundaries -currentBoundary))
    } 
  }
  mass
}


##binPeaks function! This needs the two functions above
binPeaks <- function (l, tolerance = 0.002) 
{
  nn <- sapply(l, nrow)
  nonEmpty <- nn != 0L
  samples <- rep.int(seq_along(l), nn)
  mass <- unname(unlist((lapply(l[nonEmpty], function(x) x$mz)), recursive = FALSE, use.names = FALSE))
  intensities <- unlist(lapply(l[nonEmpty], function(x) x$intensity), recursive = FALSE, use.names = FALSE)
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  mass <- binning(mass = mass, intensities = intensities, samples = samples, tolerance = tolerance)
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  lIdx <- split(seq_along(mass), samples)
  l[nonEmpty] <- mapply(FUN = function(p, i) {
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    p
  }, p = l[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  l
}

iso_hippo <- function(df, isotope_da = 1.0034, iso_tolerance = 0.0034) {
  #' iso_hippo is a function that takes 1 argument, and 2 optional ones
  #' Similarly to kopi luwak, iso hippo is a retired, senile, hippopotamus
  #' that grazes on huge chunks of data and poops isotopes and mol.ions
  #' 
  #' Ahmed's notes: Should expand on it and allow it to work with other isotopes +
  #' merge it with the repeating units function, also have it check intensity values
  #' #to make sure if its actually an isotope of the parent molecule or not
  #' 
  #' 
  #' @param df is the dataframe on which hippo eats, it has to have an mz column
  #' @param isotope_da is an optional argument, specifying the differences
  #' between mol ions and isotopes in Da
  #' @param iso_tolerance optional argument: tolerance in Da for isotope differences
  #' @importFrom tidyverse
  #' 
  #' @examples 
  #' 
  #' 
  #' @export
  mz_vector <- as.numeric(df$mz) #incase mzs are strings
  ooutput <- abs(outer(mz_vector, mz_vector, `-`))
  index_matrix <- which(ooutput <= isotope_da + iso_tolerance &
                          ooutput >= isotope_da - iso_tolerance, 
                        arr.ind = TRUE)
  iso_df <- tibble(mol_ion = mz_vector[index_matrix[, 1]],
                   isotope = mz_vector[index_matrix[, 2]],
                   diff = mol_ion - isotope)
  iso_df_f <-filter(iso_df, diff > 0)
  final_df <- select(mutate(df, isotopic_status = 
                              case_when(mz_vector %in% iso_df_f$mol_ion ~ "molecular_ion",
                                        mz_vector %in% iso_df_f$isotope ~ "c13_isotope")),
                     isotopic_status, everything())
  return(final_df)
}




# Cluster binning ---------------------------------------------------------

### Functions
# mergeData <- function(path){
#   
#   setwd(path)
#   
#   files <- list.files(path = path)
#   #files <- mixedsort(files)
#   
#   Listed <- lapply(files, function(x) read.table(x, col.names=c("mz", "intensity")))
#   
#   MyMerge <- function(x, y){
#     df <- merge(x, y, by="mz", all.x= TRUE, all.y= TRUE)
#   }
#   Listed <- Reduce(MyMerge, Listed)
#   Listed[is.na(Listed)] <- 0
#   IDvalues <- files %>% stringr::str_remove(".txt")
#   colnames(Listed) <- c("mz", IDvalues)
#   return(Listed)
# }

#binning condition. 
condition_ppm <- function (mass, tolerance) 
{
  # if the ppm between the first and the last mz of the bin is larger than the threshold ---> continue binning
  l <- length(mass)
  if (abs((mass[l]-mass[1L])/mass[l]) * 1e6 > tolerance) { 
    return(TRUE) 
  } 
  FALSE
}

condition_num <- function (mass, times) 
{
  # if the number of mz in this bin is larger than the threshold ---> continue binning
  # the floor of the times is 1 (ideally the lowest peak number in a bin = the sample number )
  # the cell of the times depends on the RAM of your computer and your expectation of running speed
  l <- length(mass)
  if (l > sample_num * times) { 
    return(TRUE) 
  } 
  FALSE
}

#Binning function. 
binning <- function (mass, tolerance, times) 
{
  n <- length(mass)
  d <- diff(mass)
  nBoundaries <- max(20L, floor(3L * log2(n))) #ensure there's enough bins
  boundary <- list(left = double(nBoundaries), right = double(nBoundaries))
  currentBoundary <- 1L
  boundary$left[currentBoundary] <- 1L
  boundary$right[currentBoundary] <- n
  
  q <- list(left = double(n), right = double(n))
  currentQ <- 1L
  
  while (currentBoundary > 0L) {
    left <- boundary$left[currentBoundary]
    right <- boundary$right[currentBoundary]
    currentBoundary <- currentBoundary - 1L
    gaps <- d[left:(right - 1L)] 
    gapIdx <- which.max(gaps) + left - 1L # find the largest gap
    
    #condition for the ending point of the left
    #two type of condition can be chosen here: based on tolerance/sample_number
    #the line below needs to be modify if you want to choose another condition
    #refer to condition_ppm() function when you need to modify
    if (condition_num(mass = mass[left:gapIdx], times)) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    }
    else {
      #mass[left:gapIdx] <- l
      q$left[currentQ] <- mass[left]
      q$right[currentQ] <- mass[gapIdx]
      currentQ <- currentQ + 1L
    }
    #condition for the ending point of the right
    if (condition_num(mass = mass[(gapIdx + 1L):right], times)) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- gapIdx + 1L
      boundary$right[currentBoundary] <- right
    }
    else {
      #mass[(gapIdx + 1L):right] <- r
      q$left[currentQ] <- mass[gapIdx]
      q$right[currentQ] <- mass[right]
      currentQ <- currentQ + 1L
    }
  }
  
  q$left <- q$left[1:currentQ - 1L]
  q$right <- q$right[1:currentQ - 1L]
  return(q)
}

ppm_calc <- function (mass1, mass2) {
  #' ppm_calc calculated the parts per million error between two different masses
  #'
  #' @examples
  #' ppm_calc(mass1, mass2)
  #'
  #' @export
  ppm_error <- ((mass1 - mass2)/mass1) * 1e6
  return(ppm_error)
}

mydensityClust <- function(data, mzMin, mzMax, dc){
  ## this is the function adapted from the densityClust() function of the 
  ## densityClust package. It returns the rho&delta decision graph and the clusters
  ## it finds. 
  ## mz need to cluster in intervals, otherwise the data amount is too large to calculate/store
  
  # pre-processing for densityClust()
  # get mz values
  mz <- data$mz
  
  # subset mz for testing
  tmz<-mz[mz>=mzMin & mz <= mzMax]
  
  # calculate the distance (in ppm) between each mz and put them in a data frame 
  dist_matrix <- as.dist(abs(usedist::dist_make(as.matrix(tmz), ppm_calc)), diag = TRUE)
  
  # distMatrix <- matrix(data = NA, nrow = length(tmz), ncol = length(tmz))
  # for (i in 1:length(tmz)){
  #   for (j in 1:length(tmz)){
  #     distMatrix[i,j] <- abs(ppm_calc(tmz[j], tmz[i]))
  #   }
  # }
  # distMatrix <- as.dist(distMatrix,diag=TRUE)
  # 
  # dc(cutoff distance): a radium within which the density will be calculated  
  tmzClust <- densityClust::densityClust(dist = dist_matrix, dc = dc, gaussian=FALSE)
  
  return(tmzClust)
}

getclusteredData <- function(data, mzClust, mzMin, mzMax){
  ## this is a function that takes the results of the densityClust() function 
  ## from densityClust package, outputs the data with rho & delta
  
  #cluster <- mzClust$clusters
  
  # find the index for sub-setting
  ind <- which(data$mz >= mzMin)
  ind1 <- ind[1]
  ind <- which(data$mz <= mzMax)
  ind2 <- length(ind)
  
  # sub-setting
  df <- data[ind1:ind2,]
  df$rho <- mzClust$rho
  df$delta <- mzClust$delta
  
  return(df)
}

clusterAssign <- function(data){
  ## this is a function takes the dataset as input, and assign clusters based 
  ## on the delta
  
  data$cluster <- rep(0, nrow(data))
  clusterID <- 1
  data$cluster[1] <- clusterID
  
  for (i in 2:nrow(data)){
    dist <- abs(data$delta[i]-data$delta[i-1])
    if(dist>5 | abs(data$rho[i] - data$rho[i-1])>2){
      clusterID <- clusterID + 1
    }
    data$cluster[i] <- clusterID
  }
  
  return (data)
} 

cluster_align <- function(bin_element){
  ## find clusters
  # rho: the minimum number of peaks within a cluster
  # delta: the minimum distance between two cluster centers
  #        (twice of the highest acceptable mass measurement error)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  
  mz_start <- bin_element[1,1]
  mz_end <- bin_element[1,2]
  
  tmz <- msdata$mz[msdata$mz >= mz_start & msdata$mz <= mz_end]
  l <- length(tmz)
  
  print(c(mz_start, l))
  
  if (l < 2){
    tmz <- mz_start
    aligned_data <- as.data.frame(t(data.frame(rep(0,sample_num+1))))
    colnames(aligned_data) <- c("aligned", file_names)
    aligned_data$aligned <- tmz
    return(aligned_data)
  }
  
  mzClust <- mydensityClust(msdata, mz_start, mz_end,  dc = 5)
  
  ## get the clustered results
  clusteredData <- getclusteredData(msdata, mzClust, mz_start, mz_end)
  
  ## store the orphaned peaks
  orphans <- clusteredData[which(clusteredData$rho == 0), ]
  
  ## find out the clusters based on delta
  # take out the orphans first
  clusteredData <- clusteredData[which(clusteredData$rho != 0), ]
  
  # cluster by the difference between the delta[i] & delta[i-1]
  if (nrow(clusteredData) != 0){
    clusteredData <- clusterAssign(clusteredData)
    
    # check the adjacent special orphan peaks
    clusterdf <- as.data.frame(table(clusteredData$cluster))
    # get the cluster ID of which cluster only contains 1 peak
    tocheck <- which(clusterdf$Freq == 1)
    # get the cluster ID of which cluster is adjacent to another cluster only contains 1 peak
    if (length(tocheck)>1){
      v <- vector()
      for (i in 2:length(tocheck)){
        dif <- tocheck[i]-tocheck[i-1]
        if (dif == 1){
          v[i] <- tocheck[i]
        }
      }
      v <- v[-which(is.na(v))]
      # get the cluster ID of which  the ppm error of these adjacent peaks is below 5
      # (check if the v is an empty vector first)
      if (is.logical(v) == FALSE ){
        v2fix <- vector()
        for (i in 1:length(v)){
          #j <- v[i]
          mz1 <- clusteredData$mz[which(clusteredData$cluster == v[i])]
          mz2 <- clusteredData$mz[which(clusteredData$cluster == v[i])-1]
          ppmerror <- ppm_calc(mz1, mz2)
          if (ppmerror < 5){
            # reassign the cluster ID for these peaks
            v2fix[i]<-v[i]
            clusteredData$cluster[which(clusteredData$cluster == v[i])] <- clusteredData$cluster[which(clusteredData$cluster == v[i])-1]
          }
        }
        v2fix <- v2fix[-which(is.na(v2fix))]
        
        # reassign the cluster ID for these peaks
        # if (is.logical(v2fix == FALSE)){
        #   for (i in 1:length(v2fix)){
        #     j <- v2fix[i]
        #     clusteredData$cluster[which(clusteredData$cluster == j)] <- j-1
        #   }
      }
    }
    
    # align the mz based on the clustering results and store them in alignedData
    alignedData <- aggregate(clusteredData$mz, by=list(cluster=clusteredData$cluster),mean)
    colnames(alignedData) <- c("cluster", "mz")
    alignedData <- alignedData[order(alignedData$mz), ]
    
    # add the aligned results to the clusteredData
    clusteredData$aligned <- rep(0, nrow(clusteredData))
    for (i in 1:nrow(clusteredData)){
      if (is.na(clusteredData$cluster[i]) == FALSE){
        clusteredData$aligned[i] <- alignedData$mz[which(clusteredData$cluster[i] == alignedData$cluster)]
      } else{
        clusteredData$aligned[i] <- clusteredData$mz[i]
      }
    }
    
  }
  
  #put the orphans back to the clusteredData
  orphans$cluster <- rep(0, nrow(orphans))
  orphans$aligned <- orphans$mz
  clusteredData <- rbind(clusteredData,orphans)
  clusteredData <- clusteredData[order(clusteredData$mz), ]
  
  # remove the rho & delta columns
  clusteredData <- subset(clusteredData, select = -c(rho, delta))
  
  ## transform the data format back to the initial 
  nulldata <- data.frame(mz = rep(0,sample_num), sample=file_names,
                         intensity = rep(0,sample_num), cluster = rep(0,sample_num), aligned = rep(0,sample_num))
  ttmpData <- rbind(nulldata, clusteredData)
  ttmpData <- ttmpData %>% pivot_wider(names_from = sample, values_from = intensity)
  ttmpData[is.na(ttmpData)] <- 0
  
  # remove the null data 
  ttmpData <- ttmpData[-1, ]
  
  # remove the sample mz & cluster ID columns
  tmpData <- subset(ttmpData, select = -c(mz, cluster))
  
  # merge the rows with same alignedmz
  # tsfData <- tmpData %>% group_by(aligned) %>% dplyr::summarize(across(.fns = max))
  
  # merge the rows with same alignedmz
  aligned_data <- tmpData %>% group_by(aligned) %>% dplyr::summarize(across(.fns = max))
  
  aligned_data
}

### Script
## merge data 
options(scipen=999)
options(digits=10)
path <- "E:/code/CTC_HU/renamed_rawdata/renamed_pos"
setwd(path)
files <- list.files(path = path)
file_names <<- files %>% stringr::str_remove(".txt")
file_list <- lapply(files, function(x) read.table(x, col.names=c("mz", "intensity")))

nn <- sapply(file_list, nrow) # each sample's row number
nonEmpty <- nn != 0L # if there's any empty sample; return logic vector
samples <- rep(file_names, nn)
mass <- unname(unlist((lapply(file_list[nonEmpty], function(x) x$mz)), recursive = FALSE, use.names = FALSE))
intensities <- unlist(lapply(file_list[nonEmpty], function(x) x$intensity), recursive = FALSE, use.names = FALSE)
s <- sort.int(mass, index.return = TRUE)
mass <- s$x
intensities <- intensities[s$ix]
samples <- samples[s$ix]

msdata <<- data.frame(mz = mass, sample = samples, intensity = intensities)


# withIDdata <- mergeData("E:/code/CTC_HU/renamed_rawdata/HB")
# 
# ## transform the data format into the one needed for densityclust()
# datalonger <- withIDdata %>% 
#   pivot_longer(!mz, names_to = "sample", values_to = "intensity")
# 
# ## remove rows with 0 intensity 
# msdata <<- datalonger[(which(datalonger$intensity != 0)), ]

## remove noise
#msdata <- noiseRemove2(msdata, estimatedNoise = 300, SNRlimit = 10)

# initialize a new data frame
sample_num <<- length(files)
# AlignedMS <- as.data.frame(t(data.frame(rep(0,sample_num+1))))
# colnames(AlignedMS) <- c("aligned", file_names)

# create bins first
bin <- binning(msdata$mz, tolerance = 2000, times = 4) #two type of conditions can be chosen; tolerance & times belong to different conditions
bin_df <- data.frame(left = bin$left, right = bin$right)
bin_df <- arrange(bin_df, bin_df$left)
bin_split <- split(bin_df, 1:nrow(bin_df))

# check the ppm between the first and last mz of the bin
# bin_df$ppm <- abs(ppm_calc(bin_df$left,bin_df$right))
# bin_df[which(bin_df$ppm<5),]

# processing without parallel
# list_files_aligned <- lapply(bin_split, cluster_align)

# processing with parallel
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl,as.list(unique(ls(environment(cluster_align)))), envir = environment(cluster_align))

list_files_aligned_pl <- parLapply(cl, bin_split, cluster_align)

stopCluster(cl)

# list to dataframe
Alligned_PM <- do.call(rbind, list_files_aligned_pl)
colnames(Alligned_PM) <- c("mz", file_names)



## this is the function adapted from the densityClust() function of the 
## densityClust package. It returns the rho&delta decision graph and the clusters
## it finds. The default rho = 0.8, delta = 0.05. These should be tuned based on
## the decision graph for each data set.

## need to cluster in intervals, otherwise the data amount is too large to calculate/store

mydensityClust <- function(data, mzMin, mzMax, dc){
  
  
  # pre-processing for densityClust()
  # get mz values
  mz <- data$mz
  
  # subset mz for testing
  tmz<-mz[mz>mzMin & mz <mzMax]
  
  # calculate the distance (in ppm) between each mz and put them in a data frame 
  distMatrix <- matrix(data = NA, nrow = length(tmz), ncol = length(tmz))
  for (i in 1:length(tmz)){
    for (j in 1:length(tmz)){
      distMatrix[i,j] <- abs(ppm_calc(tmz[j], tmz[i]))
    }
  }
  distMatrix <- as.dist(distMatrix,diag=TRUE)
  
  # inspect clustering attributes to define thresholds
  # dc(cutoff distance): a radium within which the density will be calculated  
  tmzClust <- densityClust(dist = distMatrix, dc = dc, gaussian=FALSE)
  #plot(tmzClust)
  
  
  # find clusters by assigned parameters
  #mzClust <- findClusters(tmzClust, rho=r, delta=d)
  
  return(tmzClust)
}


# Background removal_1 ----------------------------------------------------

blank_subtraction <- function(dataframe, filter_type = "median", blank_thresh = 1, nsamples_thresh = 1, blank_regx = "blank", filtered_df = FALSE){
  blanks <- select(dataframe, contains(blank_regx)) #creates a dataframe with only blanks
  samples <- select(dataframe, !contains(blank_regx)) #creates a dataframe without blanks
  
  blanks$threshold <- apply(blanks, 1, filter_type) * blank_thresh #adds column with median/max/etc. of each row 
  #times a specified threshold percentage
  
  index <- vector() #creates a vector for the for loop
  
  for(i in 1:nrow(samples)){
    samples$index[i] <- if_else(sum(samples[i, 1:ncol(samples)] <= blanks$threshold[i])/ncol(samples)*100 >= nsamples_thresh, "blank_fail", "blank_pass")
  } #calculates for each row how many data points don't exceed the initial threshold
  #if that amount meets or exceeds a certain threshold percentage of the total
  #amount of samples it will be marked as failed, if it remains below the threshold
  #percentage it will be marked as passed. Finally, a column is added to the sample 
  #dataframe that shows which m/z values passed or failed
  
  filtered_samples <- filter(samples, index == "blank_pass") #dataframe that only displays the passed values, 
  #filtering out the failed values
  print_pass <- print(nrow(filter(samples, index == "blank_pass"))) #prints the number of passed values
  print_fail <- print(nrow(filter(samples, index == "blank_fail"))) #prints the number of failed values 
  
  if(filtered_df == TRUE){
    filtered_list <- list("filtered_df" = filtered_samples, "n_passed" = print_pass, "n_failed" = print_fail)
    return(filtered_list)
  } #makes the function return a list with a dataframe containing m/z rows that 
  #passed the threshold, the number of passed values and the number of failed values
  
  if(filtered_df == FALSE){
    unfiltered_list <- list("unfiltered_df" = samples, "n_passed" = print_pass, "n_failed" = print_fail)
    return(unfiltered_list)
  } #makes the function return a list with a dataframe containing m/z rows that 
  #both passed and failed to pass the threshold, the number of passed values 
  #and the number of failed values
}



# Visualization -----------------------------------------------------------

## heatmap 
## we need transposed df   
superheat <- master_df_t %>%
  superheat(left.label.size = 0.05,
            left.label.text.size = 10,
            left.label.col = "white",
            left.label.text.angle = 90, ## the classifiers names angle
            grid.hline.col = "white",
            grid.vline.col = "white",
            grid.hline.size = 5,
            grid.vline.size = 3,
            membership.rows = classifiers,
            heat.pal = viridis::mako(100),
            pretty.order.rows = T,
            pretty.order.cols =T,
            scale = F, #scales columns
            padding = 0.1,
            legend.width=4)

# PCA and PLS-DA ----------------------------------------------------------
#visualization using factoMineR and factorextra package 
### we need transposed df   
res_pca <- PCA(master_df_pos_t, graph = FALSE) 

## Visualisation of variance explained (plot the variance against the no of dimension)
fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50),main="PCA - scree plot" )


## Extracting results of variables
var <- get_pca_var(res_pca)
fviz_pca_var(res_pca,col.var = "grey",col.circle = "grey",title="variables-PCA")

##Plotting the individuals 

fviz_pca_ind(res_pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,# Avoid text overlapping (slow if many points)
             title="individuals-PCA - names of the sample"
)

## plotting the ellipses 

fviz_pca_ind(res_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = classifiers, # color by groups
             palette = "viridis",
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Sample type",
             title = "PCA samples "
)


## different way to plot the PCA 
## PCA from basic stats 

PCA2 <- prcomp(as.matrix(master_df_pos_t), scale. = F) #PCA model using transposed df 
fviz_eig(PCA2, addlabels = TRUE, ylim = c(0, 100),main="PCA -scree plot" )
PCA_scores <- as.data.frame(PCA2$x) %>% dplyr::select(PC1, PC2)
PCA_scores$Sample <- classifiers ## we add our classifiers here 

##plotting the samples and ellipses using ggplot 
ggscatter(PCA_scores, x = "PC1", y = "PC2",
          color = "Sample", shape = "Sample", palette = "aaas",
          mean.point = TRUE, ellipse = TRUE, title = "PCA ", subtitle = "PC1(25.4%) - PC2(11.2%)") + ## add PC 1 and PC 2 percentage from scree plot
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")


#PLS-DA
## For easy visualization from mixOmics package 
## using transposed df and classifiers
plsda<-mixOmics:: plsda(master_df_pos_t, classifiers)

## plotting the samples classifiers with ellipses 
plotIndiv(plsda, ind.names = FALSE, star = TRUE, ellipse = TRUE, legend = TRUE,title = "PLS-DA samples")

## plotting the ROC curve and calculating the auc of the plsda model 
auc.plsda <- auroc(plsda,roc.comp = 1, title = "PLS-DA ROC Curve ")

## plotting the pls-da contribution to median
plotLoadings(plsda, contrib = 'max', method = 'median', comp = 1,title="pls-da contribution by median")

#volcano plot ------------------------------------------------------------

volcanoPlot <- function(data, xvalues, yvalues, title){
  ggplot(data) +
    geom_point(aes(x=xvalues, y=-log10(yvalues), colour=significant)) + ## color by significant of fdr <0.1
    ggtitle(title) +
    xlab("log2 fold change") + 
    ylab("-log10 nominal p.value") +
    #scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "right",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    scale_color_aaas() +
    theme_pubr() +
    labs_pubr()
}

## but we plot the points of log2 folchange between two groups and the values of nominal p value but it will be colored
# by significance according to adjusted p value 
## we can use here the results of wilcox or welch t test 
volcanoPlot(wilcox,wilcox$log2fc_M1_M2,wilcox$p.value, "Volcanoplo P-values M1 (welch t test and mann whiteny u) 
            significance of adjusted p values with FDR < 0.1")


# Random forest visualization  ------------------------------------------------------------------
##AUC and roc
## if we did a random forest cross validation using 5 models to obtain the final model

## and then checking the final model against a permuted model we can plot the ROC Curve of all these model to compare them against 
## each other , ideally the final model will have the highest accuracy and the permuted model will have the lowest accuracy 

roc_1<- roc(results_1$actual,results_1$prediction ) ## the actual vs the predictions of the first model .. etc 
auc_1<- auc(roc_1)

roc_2<- roc(results_2$actual,results_2$prediction )
auc_2<- auc(roc_2)

roc_3<- roc(results_3$actual,results_3$prediction )
auc_3<- auc(roc_3)

roc_4<- roc(results_4$actual,results_4$prediction )
auc_4<- auc(roc_4)

roc_5<- roc(results_5$actual,results_5$prediction )
auc_5<- auc(roc_5)

roc_final<- roc(results$actual,results$prediction )
auc_final<- auc(roc_final)

roc_permuted<- roc(results_permuted$actual,results_permuted$prediction )
auc_permuted<- auc(roc_permuted)


# ROC curve of each fold
par(mfrow=c(2,4))
plot(roc_1,col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_1, 3))))
plot(roc_2,col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_2, 3))))
plot(roc_3,col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_3, 3))))
plot(roc_4,col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_4, 3))))
plot(roc_5,col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_5, 3))))
plot(roc_final,col = "blue", main = paste("final model, AUC:", as.character(round(auc_final, 3))))
plot(roc_permuted,col = "green", main = paste("permuted model, AUC:", as.character(round(auc_permuted, 3))))


# Plot the important variables of a random forest model , importance = True in the model to plot this
varImpPlot(model,col="blue",pch= 2)

# model evaluation of random forest
## we still have to know more about this package library(rfUtilities)
## it is used from random forest cross validation 



rf_cv <- rf.crossValidation(model, training_set,seed=555, p=0.1, n=99, ntree=500) 
## use the same seed,training set , n_tree as the final model , n= number of cross validation the default is 99
#p ratio of splitting the training set for cross validation aka Proportion data withhold (default p=0.10)

# Plot cross validation versus model producers accuracy
par(mfrow=c(1,2)) 
plot(rf_cv, type = "cv", main = "CV producers accuracy")
plot(rf_cv, type = "model", main = "Model producers accuracy")

# Plot cross validation versus model oob
par(mfrow=c(1,2)) 
plot(rf_cv, type = "cv", stat = "oob", main = "CV oob error")
plot(rf_cv, type = "model", stat = "oob", main = "Model oob error")


## decision boundary from library(ElemStatLearn)
# you can do this for training set or test set .. 

set <- test_set  # or training set

## i made a model specifically for this visualization because I can't find a way around it .
# this plot uses two features only for the model or we are trying to search more about that because Ahmed wants to :D 
model_vis<-randomForest(x=cbind(training_set[,2:3]),y=training_set$samples,data=training_set,ntree=500)


X1 = seq(min(set[,2]) - 1, max(set[,2]) + 1, by = 0.01)
X2 = seq(min(set[,3]) - 1, max(set[,3]) + 1, by = 0.01)
grid_set = expand.grid(X1, X2) ## we make a grid using the sequence from min and max values of the two features we used in the model
colnames(grid_set) = colnames(set[,2:3]) ## setting the same name of the feature to our grid set
y_grid = predict(model_vis, newdata = grid_set, type = 'class') ## predict our y_grid which we will use it to color
## the background for this plot creating a decision boundary of predicition between two groups 
## using the grid set we created and our model (2 features only !)
# NOTE we need class here because we have a y_grid is a matrix!
plot(set[,2:3],
     main = 'Random Forest classification (test set)',
     xlab = 'intensity', ylab = 'intensity',
     xlim = range(X1), ylim = range(X2)) # this bit creates the limits to the values plotted this is also a part of the MAGIC as it creates the line between green and red
contour(X1, X2, matrix(as.numeric(y_grid), length(X1), length(X2)), add = TRUE)
# here we run through all the y_pred data and use if else to color the dots
# note the dots are the real data, the background is the pixel by pixel determination of prediction 
# graph the dots on top of the background give you the image
points(grid_set, pch = '.', col = ifelse(y_grid == "GroupA", 'springgreen3', 'tomato')) ## plotinng background
points(set[,2:3], pch = 21, bg = ifelse(set[, 1] == "GroupA", 'green4', 'red3')) ## plotting the real data of either training or test set



# parallel coordinates ----------------------------------------------------



## plotting the parallel coordinates from GGally package
## we use a dataset which have the group column = lipids ,, median scaled transformed intensity of each group or class as columns
ggparcoord(median_data,
           columns = 2:5, groupColumn = "Lipids",  ## columns = here we specify the range of columns of median data for each group
           showPoints = TRUE,                       ##groupColumn = either significant peaks (m/z) or lipids name
           #boxplot = TRUE,                     ## we can add boxplot to the plot
           scale="globalminmax") +              ## scale globalminmax plot the data as it is without scaling on the plot
                                                ## there is uniminmax scaling the data between 0-1
  labs(x = "groups",y= "Log2 pareto scaled intensity") +
  theme_pubr()+
  theme(panel.grid = element_blank(),panel.grid.major.x=element_line(colour="black"),
        text = element_text(size=20),  legend.position="right")+
  geom_line(size=1.5)+
  geom_point(size=4)+
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("legend.title", face = "bold",size=15)+
  font("legend.text",size=15)


## if we want to add p values to the plot and we have more than two groups so we can specify our comparison group
my_comparisons <- list( c("A", "B"), c("A", "C"), c("A", "D"),
                        c("B", "C"),c("B", "D"),c("C", "D"))
plot+stat_compare_means(comparisons = my_comparisons)+  ## we add this line to the plot
  stat_compare_means(label.y = 10)  ## computing global p value , label.y= the position of p value printed in the plot



## tagging a certain value to highlight it in the parallel coords
## here we were tagging a slope if it is positive or negative using a slope column we already had in our data
## you can tag anything in your data set
tag<-within(median_data, highlight<-if_else(slope>0, "positive", "negative"))

ggparcoord(slope_tag[order(tag$highlight),], columns = 2:5, groupColumn = "highlight", 
           ## we use the extra column we created as groupcolumn to show this highlight
           showPoints = TRUE,
           #boxplot = TRUE,
           scale="globalminmax",title="Median of logged and scaled data(pareto_scaling) - Slope highlight" ) +
  scale_color_manual(values=c("red","navy"))+
  labs(x = "groups",y= "Log2 pareto scaled intensity") +
  theme(plot.title = element_text(size=20),
        panel.grid = element_blank())+geom_line(size=1)
