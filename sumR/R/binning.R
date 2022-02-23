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
  boundary <- list(left = double(nBoundaries),
                   right = double(nBoundaries))#boundaries list
  currentBoundary <- 1L
  boundary$left[currentBoundary] <- 1L
  boundary$right[currentBoundary] <- n

  while (currentBoundary > 0L) {
    left <- boundary$left[currentBoundary]
    right <- boundary$right[currentBoundary]
    currentBoundary <- currentBoundary - 1L
    gaps <- d[left:(right - 1L)]
    ## Find the largest gap
    gapIdx <- which.max(gaps) + left - 1L

    l <- condition(mass = mass[left:gapIdx],
                   intensities = intensities[left:gapIdx],
                   samples = samples[left:gapIdx],
                   tolerance = tolerance)

    if (is.na(l[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    }
    else {
      mass[left:gapIdx] <- l
    }
    r <- condition(mass = mass[(gapIdx + 1L):right],
                   intensities = intensities[(gapIdx + 1L):right],
                   samples = samples[(gapIdx + 1L):right],
                   tolerance = tolerance)

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
      boundary$left <- c(boundary$left,
                         double(nBoundaries - currentBoundary))
      boundary$right <- c(boundary$right,
                          double(nBoundaries - currentBoundary))
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
  samples <- rep.int(seq_along(df_list),
                     sapply(df_list, nrow))

  mass <- unname(unlist((lapply(df_list[nonEmpty],
                                function(x) as.double(x$mz))),
                        recursive = FALSE, use.names = FALSE))
  intensities <- unlist(lapply(df_list[nonEmpty],
                               function(x) as.double(x$intensity)),
                        recursive = FALSE, use.names = FALSE)
  name <- unlist(lapply(df_list[nonEmpty],
                        function(x) x$name),
                 recursive = FALSE, use.names = FALSE)

  s <- sort.int(mass, index.return = TRUE) # sort vector based on masses lowest to highest
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  name <- name[s$ix]

  mass <- binning(mass = mass, intensities = intensities,
                  samples = samples, tolerance = tolerance)

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
  }, p = df_list[nonEmpty], i = lIdx, MoreArgs = NULL,
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
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
#' @param bin_plot logical value obtained from user input per default set to FALSE
#' decides if the line plot showing the reduction of the bins per iteration is shown
#' @export
iteration <- function(df_list, max_align = 8, bin_plot = F){
  count <- 0L
  bins <- c()
  df <- plyr::join_all(df_list, by = NULL,
                       type = "full", match = "all")
  df <- df[,c("name", "mz", "intensity")]
  df <- pivot_wider(df, names_from = "name",
                    id_cols = "mz",
                    values_from = "intensity")

  df <- df[order(df$mz),]
  bins <- c(bins, nrow(df))

  while (TRUE) {
    new <- binPeaks(df_list)
    count <- count + 1L
    new_df <- plyr::join_all(new, by = NULL,
                             type = "full", match = "all")
    new_df <- pivot_wider(new_df, names_from = "name",
                          id_cols = "mz",
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
  if (bin_plot == T) {
    df_bins <- as.data.frame(cbind(bins, 0:(count)))
    reduction_plot <- ggplot(df_bins) + geom_line(aes(x = V2, y = bins)) +
      labs(x = "Number of Alignments",
           y = "Number of bins")
    df <- list(df = df,
               reduction_plot = reduction_plot)
    return(df)
  }
  return(df)
}
