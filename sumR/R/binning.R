binning <- function(mass, intensities, samples, tolerance) {
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
