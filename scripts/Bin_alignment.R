#' @title binning dependency 1
#' @param mass 
#' @param intensities
#' @param samples
#' @param tolerance
condition <- function (mass, intensities, samples, tolerance){
  if (anyDuplicated(samples)){
    return(NA)
  }
  meanMass <- mean(mass)
  if (any(abs(mass - meanMass)/meanMass > tolerance)) {
    return(NA)
  }
  meanMass
}

#' @title binning dependency 2
#' @param mass 
#' @param intensities
#' @param samples
#' @param tolerance
binning <- function (mass, intensities, samples, tolerance){
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
    if (currentBoundary == nBoundaries){
      nBoundaries <- floor(nBoundaries * 1.5)
      boundary$left <- c(boundary$left, double(nBoundaries -currentBoundary))
      boundary$right <- c(boundary$right, double(nBoundaries -currentBoundary))
    } 
  }
  mass
}


##binPeaks function! This needs the two functions above
#' @title Binning of the peaks 
#' @
binPeaks <- function (l, tolerance = 0.002){
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

