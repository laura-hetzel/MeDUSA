##getting the data 

#' @title Import experimental data
#' @param files String location of the experimental files
#' @importFrom plyr 
#' @importFrom readr read_delim 
#' @export
get_data <- function(file){
  files <- list.files(path = "/Users/klarab/Documents/GitHub/sum-r/scripts", 
                      pattern = ".txt")
  for (i in seq_along(files)) {
    assign(paste("Df", i, sep = "."), read_delim(files[i], 
                                                 col_names = c("1", "mz", "2", "intensity", "3")))
  }
  l_df <- mget(ls(pattern = "Df."))
}

#' @title ppm calculation 
#' @description ppm_calc calculated the parts per million error between two different masses
#' @examples ppm_calc(mass1, mass2)
#' @export
ppm_calc <- function(mass1, mass2) {
  ppm_error <- ((mass1 - mass2)/mass1) * 1e6
  return(ppm_error)
}


#' @title alignment check 
#' @description align_check makes sure that all the m/z values are aligned/binned correctly
#' align_check takes a data frame of peaks with mz column as an argument,
#' and the coordinates for the plot to be zoomed in on, as an optional argument
#' align_check outputs a list of three elements:
#' 1- Dataframe of 1 column containing the ppm error values
#' 2- boxplot of the ppm erro values with xcoords zoomed in to -20,0 (default)
#' 3- table of summary stats of ppm error values
#' @importFrom dplyr
#' @importFrom tibble as_tibble_col
#' @importFrom ggplot2 ggplot
#' @export
align_check <- function(data_frame_fn, xcoords = c(-50, 0)) {
  odd_ind_fn <- seq(3, length(data_frame_fn$mz), 2)
  even_ind_fn <- seq(2, length(data_frame_fn$mz),2)
  ppm_err_fn <- ppm_calc(data_frame_fn$mz[even_ind_fn], data_frame_fn$mz[odd_ind_fn]) %>% 
    as_tibble_col(column_name = "ppm_error") 
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


#' @title Deletion of unwanted samples  
#' @description binning dependency 1
#' delete any duplication of a samples 
#' delete any sample that is kicked out by the tolerance 
#' @param mass 
#' @param samples
#' @param intensities
#' @param tolerance
condition <- function(mass, intensities, samples, tolerance) {
  if (anyDuplicated(samples)) {
    return(NA)
  }
  if (any(abs(mass - mean(mass))/mean(mass) > tolerance)) {
    return(NA)
  }
  return(mean(mass))
}

#' @title Setting of the bin boundaries 
#' @description binning dependency 2
#' @param mass 
#' @param intensities
#' @param samples
#' @param tolerance
binning <- function(mass, intensities, samples, tolerance){
  n <- length(mass) #number of peaks 
  d <- diff(mass)   #difference between masses 
  nBoundaries <- max(20L, floor(3L * log(n)))#setting of amount of bins
  boundary <- list(left = double(nBoundaries), right = double(nBoundaries))#boundaries list 
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
      boundary$left <- c(boundary$left, double(nBoundaries - currentBoundary))
      boundary$right <- c(boundary$right, double(nBoundaries - currentBoundary))
    } 
  }
  return(mass)
}

#' @title Binning of the peaks 
#' @description binPeaks function! This needs the two functions above
#' @param list of dataframes 
binPeaks <- function(l, tolerance = 5e-6) {
  nonEmpty <- sapply(l, nrow) != 0L #checking if the list is not empty
  samples <- rep.int(seq_along(l), sapply(l, nrow))
   mass <- unname(unlist((lapply(l[nonEmpty], function(x) as.double(x$mz))), recursive = FALSE, use.names = FALSE))
   intensities <- unlist(lapply(l[nonEmpty], function(x) as.double(x$intensity)), recursive = FALSE, use.names = FALSE)
   s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest 
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  mass <- binning(mass = mass, intensities = intensities, samples = samples, tolerance = tolerance)
  s <- sort.int(mass, index.return = TRUE)# sort results into mass lowest to highest 
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  lIdx <- split(seq_along(mass), samples)
  l[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses 
    p = NULL
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    return(as.data.frame(p))
  }, p = l[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  return(l)
}

#' @title Merge dataframes of all aligned samples 
#' @param dataframe obtained from experimental data using `binPeaks`
#' @param path obtained from user input -> path to where the files are stored
createDF <- function(l) {
  MyMerge <- function(x, y){
    df <- merge(x, y, by = "mz", all.x = TRUE, all.y = TRUE)
  }
  
  l <- Reduce(MyMerge, l)
  l[is.na(l)] <- 0
  IDvalues <- pm_files %>% stringr::str_remove(".txt")
  colnames(l) <- c("mz", IDvalues)
  return(l)
}

#-------------------------------------------
##testing  
#-------------------------------------------
l <- get_data(file)
pm_files <- list.files(path = "/Users/klarab/Documents/GitHub/sum-r/scripts", pattern = ".txt")

for (i in 1:2) {
  l <- binPeaks(l)
}
count <- 2L
while (TRUE) {
  new <- binPeaks(l)
  count <- count + 1L
  new_df <- createDF(new)
  df <- createDF(l)
  if (nrow(df) == nrow(new_df)) {
    break
  }
  if (count == 10) {
    break
  }
  l <- new
}

align_check(l)



write.table(t, file = "0_iteration.txt", sep = "\t", row.names = FALSE , col.names = F)

getAnywhere("MALDIQuant")



