

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
align_check <- function(data_frame_fn, xcoords = c(-20, 0)) {
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
      boundary$left <- c(boundary$left, double(nBoundaries - currentBoundary))
      boundary$right <- c(boundary$right, double(nBoundaries - currentBoundary))
    } 
  }
  return(mass)
}

#' @title Binning of the peaks 
#' @description binPeaks function! This needs the two functions above
#' @param list of dataframes 
binPeaks <- function(l, ppm = 5) {
  nn <- sapply(l, nrow)
  nonEmpty <- nn != 0L
  samples <- rep.int(seq_along(l), nn)
   mass <- unname(unlist((lapply(l[nonEmpty], function(x) x$mz)), recursive = FALSE, use.names = FALSE))
  intensities <- unlist(lapply(l[nonEmpty], function(x) x$intensity), recursive = FALSE, use.names = FALSE)
  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest 
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  mass <- binning(mass = mass, intensities = intensities, samples = samples, tolerance = ppm)
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  lIdx <- split(seq_along(mass), samples)
  l[nonEmpty] <- mapply(FUN = function(p, i) {
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    return(p)
  }, p = l[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  return(l)
}

#-------------------------------------------
##testing 
#-------------------------------------------

binPeaks(l)

t <- binPeaks(l)
t <- l$df1

align_check(t)




##getting the data 

#' @title Import experimental data
#' @param file String location of the experimental file
#' @importFrom stats na.omit
#' @importFrom readr read_delim
#' @export
get_data <- function(file){
  # Downloading the experimental data and converting it to a data frame
  ex_data_df <- read_delim(file, col_names = T)
  colnames(ex_data_df) <- c("x","mz","y", "intensity", "z")
  ex_data_df <- ex_data_df[,c("mz","intensity")]
  return(ex_data_df)
}

create_df_list <- function(df1,df2, df3, df4, df5, df6, df7, df8, df9, df10){
  l <- mget(ls(pattern = "df[1,2,3,4,5,6,7,8,9,10]"))
}
  

dfs <- function(file) {
  df1 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test1.txt")
  df2 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test2.txt")
  df3 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test3.txt")
  df4 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test4.txt")
  df5 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test5.txt")
  df6 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test6.txt")
  df7 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test7.txt")
  df8 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test8.txt")
  df9 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test9.txt")
  df10 <- get_data("/Users/klarab/Documents/GitHub/sum-r/scripts/bin_test/test10.txt")
  l <- create_df_list(df1,df2, df3, df4, df5, df6, df7, df8, df9, df10)
}
l <- dfs(file)




