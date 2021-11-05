#' ppm_calc calculated the parts per million error between two different masses
#'
#' @examples
#' ppm_calc(mass1, mass2)
#'
#' @export
ppm_calc <- function (mass1, mass2) {
  ppm_error <- ((mass1 - mass2)/mass1) * 1e6
  return(ppm_error)
}


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
#' @import ppm_calc
#' @examples
#'
#'
#' @export
align_check <- function(data_frame_fn, xcoords = c(-20, 0)) {
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

condition <- function (mass, intensities, samples, tolerance) {
  if (anyDuplicated(samples)) {
    return(NA)
  }
  meanMass <- mean(mass)
  if (any(abs(mass - meanMass)/meanMass > tolerance)) {
    return(NA)
  }
  meanMass
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
