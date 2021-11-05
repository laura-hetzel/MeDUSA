#' @import usedist
#' @import densityClust
mydensityClust_1 <- function(data, mzMin, mzMax, dc) {
  ## this is the function adapted from the densityClust() function of the
  ## densityClust package. It returns the rho&delta decision graph and the clusters
  ## it finds.
  ## mz need to cluster in intervals, otherwise the data amount is too large to calculate/store

  # pre-processing for densityClust()

  # subset mz for testing
  tmz <- dplyr::filter(data, mz >= mzMin & mz <= mzMax) %>% pull(mz)

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
  densityClust::densityClust(dist = dist_matrix,
                             dc = dc, gaussian = FALSE)
}

getclusteredData <- function(data, mzClust, mzMin, mzMax){
  ## this is a function that takes the results of the densityClust() function
  ## from densityClust package, outputs the data with rho & delta

  #cluster <- mzClust$clusters

  # find the index for sub-setting
  # ind <- which(data$mz >= mzMin)
  # ind1 <- ind[1]
  # ind <- which(data$mz <= mzMax)
  # ind2 <- length(ind)

  df <- data %>%
    filter(mz >= mzMin & mz <= mzMax) %>%
    mutate(rho = mzClust$rho,
           delte = mzClust$delta)

  # # sub-setting
  # df <- data[ind1:ind2,]
  # df$rho <- mzClust$rho
  # df$delta <- mzClust$delta

  return(df)
}

clusterAssign <- function(data){
  ## this is a function takes the dataset as input, and assign clusters based
  ## on the delta

  data$cluster <- rep(0, nrow(data))
  clusterID <- 1
  data$cluster[1] <- clusterID

  for (i in 2:nrow(data)) {
    dist <- abs(data$delta[i] - data$delta[i - 1])
    if (dist > 5 | abs(data$rho[i] - data$rho[i - 1]) > 2) {
      clusterID <- clusterID + 1
    }
    data$cluster[i] <- clusterID
  }

  return(data)
}

cluster_align <- function(bin_element){
  ## find clusters
  # rho: the minimum number of peaks within a cluster
  # delta: the minimum distance between two cluster centers
  #        (twice of the highest acceptable mass measurement error)

  mz_start <- bin_element[1, 1]
  mz_end <- bin_element[1, 2]

  tmz <- msdata$mz[msdata$mz >= mz_start & msdata$mz <= mz_end]
  l <- length(tmz)

  print(c(mz_start, l))

  if (l < 2) {
    tmz <- mz_start
    aligned_data <- as.data.frame(t(data.frame(rep(0, sample_num + 1))))
    colnames(aligned_data) <- c("aligned", file_names)
    aligned_data$aligned <- tmz
    return(aligned_data)
  }

  mzClust <- mydensityClust(msdata, mz_start, mz_end,  dc = 5) #_1 or _2 ???

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
    if (length(tocheck) > 1) {
      v <- vector()
      for (i in 2:length(tocheck)) {
        dif <- tocheck[i] - tocheck[i - 1]
        if (dif == 1) {
          v[i] <- tocheck[i]
        }
      }
      v <- v[-which(is.na(v))]
      # get the cluster ID of which  the ppm error of these adjacent peaks is below 5
      # (check if the v is an empty vector first)
      if (is.logical(v) == FALSE) {
        v2fix <- vector()
        for (i in 1:length(v)) {
          #j <- v[i]
          indexes <- which(clusteredData$cluster == v[i])
          mz1 <- clusteredData$mz[indexes]
          mz2 <- clusteredData$mz[indexes - 1]
          ppmerror <- ppm_calc(mz1, mz2)
          if (ppmerror < 5) {
            # reassign the cluster ID for these peaks
            v2fix[i] <- v[i]
            clusteredData$cluster[indexes] <- clusteredData$cluster[indexes - 1]
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
    alignedData <- aggregate(clusteredData$mz,
                             by = list(cluster = clusteredData$cluster),
                             mean)
    colnames(alignedData) <- c("cluster", "mz")
    alignedData <- alignedData[order(alignedData$mz), ]

    # add the aligned results to the clusteredData
    clusteredData$aligned <- rep(0, nrow(clusteredData))
    for (i in 1:nrow(clusteredData)) {
      if (!is.na(clusteredData$cluster[i])) {
        indexes <- which(clusteredData$cluster[i] == alignedData$cluster)
        clusteredData$aligned[i] <- alignedData$mz[indexes]
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
  nulldata <- data.frame(mz = rep(0, sample_num),
                         sample = file_names,
                         intensity = rep(0, sample_num),
                         cluster = rep(0, sample_num),
                         aligned = rep(0, sample_num))
  ttmpData <- rbind(nulldata, clusteredData)
  ttmpData <- ttmpData %>% tidyr::pivot_wider(names_from = sample,
                                              values_from = intensity)
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


## this is the function adapted from the densityClust() function of the
## densityClust package. It returns the rho&delta decision graph and the clusters
## it finds. The default rho = 0.8, delta = 0.05. These should be tuned based on
## the decision graph for each data set.

## need to cluster in intervals, otherwise the data amount is too large to calculate/store

mydensityClust_2 <- function(data, mzMin, mzMax, dc) {

  # pre-processing for densityClust()
  # get mz values
  mz <- data$mz

  # subset mz for testing
  tmz <- mz[mz > mzMin & mz < mzMax]

  # calculate the distance (in ppm) between each mz and put them in a data frame
  distMatrix <- matrix(data = NA, nrow = length(tmz), ncol = length(tmz))
  for (i in 1:length(tmz)) {
    for (j in 1:length(tmz)) {
      distMatrix[i,j] <- abs(ppm_calc(tmz[j], tmz[i]))
    }
  }
  distMatrix <- as.dist(distMatrix,diag = TRUE)

  # inspect clustering attributes to define thresholds
  # dc(cutoff distance): a radium within which the density will be calculated
  tmzClust <- densityClust(dist = distMatrix, dc = dc, gaussian = FALSE)
  #plot(tmzClust)


  # find clusters by assigned parameters
  #mzClust <- findClusters(tmzClust, rho=r, delta=d)

  return(tmzClust)
}


