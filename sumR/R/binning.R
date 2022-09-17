#' @title Deletion of unwanted samples
#' @description binning dependency 1
#' delete any duplication of a samples
#' delete any sample that is kicked out by the tolerance
#' @param mass obtained from the input files
#' @param intensities obtained from the input files
#' @param samples obtained from the input files
#' @param tolerance (optional) obtained from the user input,
#' default value 5e-6
condition <- function(mass, intensities, samples, tolerance = 5e-6) {
  if (anyDuplicated(samples)) {
    return(NA)
  }
  m <- mean(mass)
  if (any(abs(mass - m) / m > tolerance)) {
    return(NA)
  }
  return(m)
}

#' @title Setting of the bin boundaries
#' @description binning dependency 2
#' @param mass obtained from the input files
#' @param intensities obtained from the input files
#' @param samples obtained from the input files
#' @param tolerance (optional) obtained from the user input,
#' default value set to 5e-6
binning <- function(mass, intensities, samples, tolerance = 5e-6) {
  n <- length(mass) # number of peaks
  d <- diff(mass) # difference between masses
  nBoundaries <- max(20L, floor(3L * log(n))) # setting of amount of bins
  boundary <- list(
    left = double(nBoundaries),
    right = double(nBoundaries)
  ) # boundaries list
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

    l <- condition(
      mass = mass[left:gapIdx],
      intensities = intensities[left:gapIdx],
      samples = samples[left:gapIdx],
      tolerance = tolerance
    )

    if (is.na(l[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- left
      boundary$right[currentBoundary] <- gapIdx
    } else {
      mass[left:gapIdx] <- l
    }
    r <- condition(
      mass = mass[(gapIdx + 1L):right],
      intensities = intensities[(gapIdx + 1L):right],
      samples = samples[(gapIdx + 1L):right],
      tolerance = tolerance
    )

    if (is.na(r[1L])) {
      currentBoundary <- currentBoundary + 1L
      boundary$left[currentBoundary] <- gapIdx + 1L
      boundary$right[currentBoundary] <- right
    } else {
      mass[(gapIdx + 1L):right] <- r
    }
    if (currentBoundary == nBoundaries) {
      nBoundaries <- floor(nBoundaries * 1.5)
      boundary$left <- c(
        boundary$left,
        double(nBoundaries - currentBoundary)
      )
      boundary$right <- c(
        boundary$right,
        double(nBoundaries - currentBoundary)
      )
    }
  }
  return(mass)
}
