#' @title Import experimental data
#' @param files String location of the experimental files
#' @importFrom readr read_delim 
#' @export
get_data <- function(files){
  files <- list.files(path = "/Users/klarab/Documents/GitHub/sum-r/scripts", 
                      pattern = ".txt")
  for (i in seq_along(files)) {
    assign(paste("Df", i, sep = "."), 
           read_delim(files[i],
                      col_types = cols(mz = col_double(), 
                                      intensity = col_double()), 
                      id = "name", 
                      col_names = (c("1", "mz", "2", "intensity", "3"))))
                                                
    
  }
  df_list <- mget(ls(pattern = "Df."))
  return(df_list)
}

#' @title ppm calculation 
#' @description ppm_calc calculated the parts per million error between two different masses
#' @param mass1 obtained from the input files
#' @param mass2 obtained from the input files
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
#' @importFrom dplyr %>%
#' @importFrom tibble as_tibble_col
#' @param aligned_df dataframe of aligned peaks obtained iteration function 
#' @param boxplot logical value deciding output of ppm error boxplot per default set to FALSE 
#' @export
align_check <- function(data_frame_fn) {
  odd_ind_fn <- seq(3, length(data_frame_fn$mz), 2)
  even_ind_fn <- seq(2, length(data_frame_fn$mz),2)
  ppm_err_fn <- ppm_calc(data_frame_fn$mz[even_ind_fn], data_frame_fn$mz[odd_ind_fn]) %>% 
    as_tibble_col(column_name = "ppm_error") #just select the column in the df 
    #ppm_err_fn$ppm_error 
  ppm_err_summary_fn <- summary(ppm_err_fn)
  ppm_err_list <-  list(ppm_error_df = ppm_err_fn,
                        summary_stats =  ppm_err_summary_fn
  )
  return(ppm_err_list)
}

#' @title Boxplot of the ppm errors 
#' @param ppm_err_fn dataframe obtained from `align_check`
#' @param xcoords vector obtained from user input or use of default value c(-50, 0)
#' @importFrom ggplot2 ggplot
ppm_err_plot <- function(ppm_err_fn, xcoords = c(-50, 0)){
  ppm_err_plot_fn <- ggplot(ppm_err_fn, aes(x = ppm_error)) +
    geom_boxplot() +
    coord_cartesian(xlim = xcoords) +
    ggtitle("ppm error boxplot") +
    theme_classic(base_size = 20)
  return(ppm_err_plot_fn)
}

#' @title Checking the results of the alignment with boxplot output if desired
#' @param aligned_peaks dataframe obtained from `iteration`
#' @param boxplot logical value obtained from user input per default set to FALSE
#' @importFrom ggplot2 ggplot
check_process <- function(aligned_peaks, boxplot = F){
  check <- align_check(aligned_peaks)
  if (boxplot == T) {
    error_plot <- ppm_err_plot(check[["ppm_error_df"]])
    check[[3]] <- error_plot
  }
  return(check)
}

#' @title Deletion of unwanted samples  
#' @description binning dependency 1
#' delete any duplication of a samples 
#' delete any sample that is kicked out by the tolerance 
#' @param mass obtained from the input files 
#' @param samples obtained from the input files
#' @param intensities obtained from the input files
#' @param tolerance obtained from the user input or use of default value 5e-6
condition <- function(mass, intensities, samples, tolerance = 5e-6) {
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
#' @param mass obtained from the input files
#' @param intensities obtained from the input files
#' @param samples obtained from the input files
#' @param tolerance obtained from the user input or use of default value 5e-6
binning <- function(mass, intensities, samples, tolerance = 5e-6){
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
#' @param df_list list of dataframes obtained from the input files 
#' @param tolerance obtained from the user input or use of default value 5e-6
binPeaks <- function(df_list, tolerance = 5e-6) {
  nonEmpty <- sapply(df_list, nrow) != 0L #checking if the list is not empty
  samples <- rep.int(seq_along(df_list), sapply(df_list, nrow))
   mass <- unname(unlist((lapply(df_list[nonEmpty], function(x) as.double(x$mz))), recursive = FALSE, use.names = FALSE))
   intensities <- unlist(lapply(df_list[nonEmpty], function(x) as.double(x$intensity)), recursive = FALSE, use.names = FALSE)
   name <- unlist(lapply(df_list[nonEmpty], function(x) x$name), recursive = FALSE, use.names = FALSE)
   s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest 
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  name <- name[s$ix]
  mass <- binning(mass = mass, intensities = intensities, samples = samples, tolerance = tolerance)
  s <- sort.int(mass, index.return = TRUE)# sort results into mass lowest to highest 
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  name <- name[s$ix]
  lIdx <- split(seq_along(mass), samples)
  df_list[nonEmpty] <- mapply(FUN = function(p, i) { # reassigning of the new masses 
    p = NULL
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    p$name <- name[i]
    return(as.data.frame(p))
  }, p = df_list[nonEmpty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  return(as.list(df_list))
}

#' @title iterations of the alignment function over the data 
#' @description the code iterates over the data as long as 
#' the number of bins changes, when the number of bins 
#' doesn't change anymore iterations stop, maximal number
#' of iterations is set to 15
#' Also, plotting the decrease of the bins 
#' @importFrom plyr join_all
#' @importFrom tidyr pivot_wider
#' @param df_list list of dataframes obtained from `binPeaks` 
#' @param max_align value obtained from user input or use of default value 8
#' @export
iteration <- function(df_list, max_align = 8){
  count <- 0L
  bins <- c()
  df <- plyr::join_all(df_list, by = NULL, type = "full", match = "all")
  df <- df[,c("name", "mz", "intensity")]
  df <- pivot_wider(df, names_from = "name", id_cols = "mz", 
                    values_from = "intensity")
  df <- df[order(df$mz),]
  bins <- c(bins, nrow(df))
  while (TRUE) {
    new <- binPeaks(df_list)
    count <- count + 1L
    new_df <- plyr::join_all(new, by = NULL, type = "full", match = "all")
    new_df <- pivot_wider(new_df, names_from = "name", id_cols = "mz", 
                          values_from = "intensity")
    new_df <- new_df[order(new_df$mz),]
    bins <- c(bins, nrow(new_df))
    if (nrow(df) == nrow(new_df)) {
      break
    }
    if (count == max_align) {
      break
    }
    df_list <- new
    df <- new_df
  }
  return(df)
}

#' @title Plotting of the decrease in peaks 
#' @param df_bins dataframe obtained from `iteration`
#' @param bins vector containing number of bins obtained from `iteration`
#' @param count value representing the number of iterations obtained from `iteration`
#' @importFrom ggplot2 ggplot
bin_plot <- function(df_bins, bins, count){
  df_bins <- as.data.frame(cbind(bins, 0:(count)))
  reduction_plot <- ggplot(df_bins) + geom_line(aes(x = V2, y = bins)) +
    labs(x = "Number of Alignments",
         y = "Number of bins") 
  return(reduction_plot)
}

#-------------------------------------------
##testing  
#-------------------------------------------
df_list <- get_data(file)
aligned_peaks <- iteration(df_list, max_align = 5)
align_analysis <- check_process(aligned_peaks, boxplot = T)
