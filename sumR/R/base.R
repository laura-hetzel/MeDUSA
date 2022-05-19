#' @title (LEGACY, use read_mzml instead) Read MS data
#' @param path path to the file
#' @return Data.frame with ordered mz values
#' @export
#' @importFrom stringr str_remove
#' @importFrom dplyr bind_rows
#' @importFrom utils read.table
#' @importFrom stats setNames
read_msdata <- function(path = "data") {
  files <- list.files(path, full.names = T)
  file_names <- str_remove(string = files, pattern = ".txt")
  file_list <- lapply(setNames(files, file_names), function(x) {
    read.table(x, col.names = c("mz", "intensity"))
  })
  nn <- sapply(file_list, nrow) # each sample's row number
  non_empty <- nn != 0L # if there's any empty sample; return logic vector
  ms_data <- bind_rows(file_list[non_empty], .id = "sample")

  return(ms_data[order(ms_data$mz), ])
}

#' @title Convert RAW files to mzML
#' @param folder
#' @param outputFolder
#' @param options
#' @export
rawToMzml <- function(folder, output = getwd(), rt = NULL, options = ""){

  options <- paste(options, collapse = " ")
  files <- list.files(file.path(folder), full.names = T, pattern = "^(?!.*mzML).*")
  output <- file.path(output)
  if (!dir.exists(output)) dir.create(output)
  key <- "Directory\\shell\\Open with MSConvertGUI\\command"
  reg <- tryCatch(utils::readRegistry(key, "HCR"), error = function(e) NULL)
  msconvert <- file.path(dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])), "msconvert.exe")
  if (is.null(msconvert)) return(NULL)

  rt <- ifelse(is.null(rt), "", sprintf('--filter "scanTime [%s, %s]"', rt[1], rt[2]))

  pbapply::pblapply(files, function(file){
    command <- sprintf('"%s" "%s" %s %s -o "%s"',
                       file.path(msconvert), file, rt, options, output)
    system(command, show.output.on.console = F)
  })
  list.files(folder, full.names = T, pattern = ".mzML")
}


#' @title ppm calculation
#' @description ppm_calc calculated the parts per million error between two different masses
#' @param mass1 input from the `align_check` functions
#' @param mass2 input from the `align_check` functions
ppm_calc <- function(mass1, mass2) {
  ppm_error <- ((mass1 - mass2) / mass1) * 1e6
  return(ppm_error)
}


#' @title alignment check
#' @param aligned_peaks dataframe of aligned peaks obtained iteration function
#' @export
align_check <- function(aligned_peaks) {
  odd_ind_fn <- seq(3, length(aligned_peaks$mz), 2)
  even_ind_fn <- seq(2, length(aligned_peaks$mz), 2)
  ppm_err_fn <- data.frame("ppm_error" = ppm_calc(aligned_peaks$mz[even_ind_fn], aligned_peaks$mz[odd_ind_fn]))
  return(ppm_err_fn)
}

#' @title Boxplot of the ppm errors
#' @param ppm_err_fn dataframe obtained from `align_check`
#' @importFrom ggplot2 ggplot
#' @export
ppm_err_plot <- function(ppm_err_fn) {
  ppm_err_plot_fn <- ggplot(ppm_err_fn, aes(x = ppm_error)) +
    geom_boxplot() +
    ggtitle("ppm error boxplot") +
    theme_classic(base_size = 20)
  return(ppm_err_plot_fn)
}

#' @title Checking the results of the alignment with boxplot output if desired
#' @description here the user can choose what kind of analysis they want to have on
#' their alignment, check_process makes sure that all the m/z values are aligned/binned
#' correctly, check_process takes a data frame of peaks with mz column as an argument,
#' and the coordinates for the plot to be zoomed in on, as an optional argument
#' check_process outputs either a dataframe (1) or a list of two to three elements:
#' 1- Dataframe of 1 column containing the ppm error values
#' 2- (optional)table of summary stats of ppm error values
#' 3- (optional)- boxplot of the ppm erro values with xcoords zoomed in to -50,0 (default)
#' @param aligned_peaks dataframe obtained from `iteration`
#' @param summary_errors logical value obtained from user input per default set to FALSE
#' @param boxplot logical value obtained from user input per default set to FALSE
#' @param xcoords vector obtained from user input or use of default value c(-50, 0)
#' @importFrom ggplot2 ggplot
#' @export
bin_align_check_process <- function(aligned_peaks, summary_errors = F,
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


#' @title Deletion of unwanted samples
#' @description binning dependency 1
#' delete any duplication of a samples
#' delete any sample that is kicked out by the tolerance
#' @param mass obtained from the input files
#' @param intensities obtained from the input files
#' @param samples obtained from the input files
#' @param tolerance obtained from the user input or use of default value 5e-6
condition <- function(mass, intensities, samples, tolerance = 5e-6) {
  if (anyDuplicated(samples)) {
    return(NA)
  }
  if (any(abs(mass - mean(mass)) / mean(mass) > tolerance)) {
    return(NA)
  }
  return(mean(mass))
}


#' @title Background removal
#' @param dataframe
#' @param filter_type
#' @param blank_thresh
#' @param nsamples_thresh
#' @param blank_regx
#' @param filtered_df
#' @importFrom dplyr select contains if_else
#' @export
blank_subtraction <- function(dataframe, filter_type = "median",
                              blank_thresh = 1, nsamples_thresh = 1,
                              blank_regx = "blank", filtered_df = FALSE) {
  ## Creates a dataframe with only blanks
  blanks <- select(dataframe, contains(blank_regx))
  ## Creates a dataframe without blanks
  samples <- select(dataframe, !contains(blank_regx))

  ## Adds column with median/max/etc. of each row times a specified
  ## threshold percentage
  blanks$threshold <- apply(blanks, 1, filter_type) * blank_thresh

  ## Calculates for each row how many data points don't exceed the initial
  ## threshold if that amount meets or exceeds a certain threshold percentage
  ## of the total amount of samples it will be marked as failed, if it remains
  ## below the threshold percentage it will be marked as passed. Finally, a
  ## column is added to the sample dataframe that shows which m/z values passed
  ## or failed
  for (i in seq_len(nrow(samples))) {
    below_thresh <- sum(samples[i, seq_len(ncol(samples))] <=
      blanks$threshold[i])
    boundary <- below_thresh / ncol(samples) * 100 >= nsamples_thresh
    samples$index[i] <- if_else(boundary, "blank_fail", "blank_pass")
  }

  ## Dataframe that only displays the passed values, filtering out failed values

  filtered_samples <- samples[samples$index == "blank_pass", ]
  ## Prints the number of passed values
  print_pass <- sum(samples == "blank_pass")
  ## Prints the number of failed values
  print_fail <- sum(samples == "blank_fail")

  if (filtered_df == TRUE) {
    ## Makes the function return a list with a dataframe containing m/z
    ## rows that passed the threshold, the number of passed values and the
    ## number of failed values
    filtered_list <- list(
      "filtered_df" = filtered_samples,
      "n_passed" = print_pass,
      "n_failed" = print_fail
    )
    return(filtered_list)
  }

  ## Makes the function return a list with a dataframe containing m/z rows that
  ## both passed and failed to pass the threshold, the number of passed values
  ## and the number of failed values
  unfiltered_list <- list(
    "unfiltered_df" = samples,
    "n_passed" = print_pass,
    "n_failed" = print_fail
  )
  return(unfiltered_list)
}

#' @title Calculate nominal mass from a formula
#' @param formulas formulas in string format like 'c6h12o6'
#' @importFrom stringr str_extract_all
#' @importFrom PeriodicTable mass
calculate_nominal_mass <- function(formulas) {
  sapply(formulas, function(formula) {
    vec <- as.vector(str_extract_all(formula, "([A-Za-z][[:digit:]]+|[A-Z][a-z]{0,1})",
      simplify = T
    ))
    sum(sapply(vec, function(comb) {
      if (!grepl("[[:digit:]]", comb)) {
        atom <- comb
        mult <- 1
      } else {
        atom <- str_extract(comb, "[[:alpha:]]*")
        mult <- as.integer(str_extract(comb, "[[:digit:]].*"))
      }

      sum(round(PeriodicTable::mass(atom)) * mult)
    }))
  })
}
