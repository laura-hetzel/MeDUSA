#' @title ppm calculation
#' @description ppm_calc calculates the parts per million
#' error between two adjacent masses
#' @param mass1 input from the `align_check` functions
#' @param mass2 input from the `align_check` functions
ppm_calc <- function(mass1, mass2) {
  ppm_error <- ((mass1 - mass2) / mass1) * 1e6
  return(ppm_error)
}


#' @title alignment quality control
#' @description uses the ppm_calc function to check
#' ppm error of aligned peaks
#' @param aligned_peaks dataframe of aligned peaks
#' obtained from `iteration`
#' @export
align_check <- function(aligned_peaks) {
  odd_ind_fn <- seq(3, length(aligned_peaks$mz), 2)
  even_ind_fn <- seq(2, length(aligned_peaks$mz), 2)
  ppm_err_fn <- data.frame("ppm_error" = ppm_calc(aligned_peaks$mz[even_ind_fn], aligned_peaks$mz[odd_ind_fn]))
  return(ppm_err_fn)
}

#' @title Boxplot of the ppm errors
#' @param ppm_err_fn dataframe obtained from `align_check`
#' @importFrom ggplot2 ggplot .data geom_boxplot ggtitle theme_classic
ppm_err_plot <- function(ppm_err_fn) {
  ggplot(ppm_err_fn, aes(x = .data$ppm_error)) +
    geom_boxplot() +
    ggtitle("ppm error boxplot") +
    theme_classic(base_size = 20)
}

#' @title Alignment analysis (optional summary & boxplot of ppm errors)
#' @description various kinds of visualizations are possible,
#' check_binning ensures m/z values are correctly aligned/binned,
#' check_binning outputs either a dataframe (1) or
#' a list of two to three elements:
#' 1- Dataframe of 1 column containing the ppm error values
#' 2- (optional)table of summary stats of ppm error values
#' 3- (optional)boxplot of the ppm error values with xcoords
#' zoomed in to -50,0 (default)
#' decide via logical (optional) values which results are
#' created, if there is no input only 1 will be created
#' @param aligned_peaks dataframe obtained from `iteration`
#' @param summary_errors (optional)logical value obtained from
#' user input, per default set to FALSE
#' @param boxplot (optional) logical value obtained from user
#' input, per default set to FALSE
#' @param xcoords (optional) vector defining zoom on boxplot
#' obtained from user input, default value set to c(-50, 0)
#' @importFrom ggplot2 ggplot
#' @export
check_binning <- function(aligned_peaks, summary_errors = F,
                                    boxplot = F, xcoords = c(-50, 0)) {
  check <- align_check(aligned_peaks)
  x <- 2
  if (summary_errors | boxplot == T) {
    input <- check
    check <- as.list(check)
  }
  if (summary_errors == T) {
    ppm_err_summary_fn <- summary(input)
    check[[x]] <- ppm_err_summary_fn
    x <- x + 1
  }
  if (boxplot == T) {
    error_plot <- ppm_err_plot(input) +
      coord_cartesian(xlim = xcoords)
    check[[x]] <- error_plot
  }
  return(check)
}







