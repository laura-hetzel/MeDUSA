#' Binning function
#'
#' @param mass
#' @param tolerance
#' @param times
#' @param sample_num
#'
#' @return
#' @export
#'
#' @examples
binning <- function(mass, tolerance, times, sample_num) {
  n <- length(mass)
  d <- diff(mass)
  n_boundaries <- max(20L, floor(3L * log2(n))) #ensure there's enough bins
  boundary <- list(left = double(n_boundaries), right = double(n_boundaries))
  current_boundary <- 1L
  boundary$left[current_boundary] <- 1L
  boundary$right[current_boundary] <- n

  q <- list(left = double(n), right = double(n))
  current_q <- 1L

  while (current_boundary > 0L) {
    left <- boundary$left[current_boundary]
    right <- boundary$right[current_boundary]
    current_boundary <- current_boundary - 1L
    gaps <- d[left:(right - 1L)]
    ## Find the largest gap
    gap_idx <- which.max(gaps) + left - 1L

    #condition for the ending point of the left
    #two type of condition can be chosen here: based on tolerance/sample_number
    #the line below needs to be modify if you want to choose another condition
    #refer to condition_ppm() function when you need to modify
    if (condition_num(mass = mass[left:gap_idx],
                      times, sample_num)) {
      current_boundary <- current_boundary + 1L
      boundary$left[current_boundary] <- left
      boundary$right[current_boundary] <- gap_idx
    }
    else {
      #mass[left:gap_idx] <- l
      q$left[current_q] <- mass[left]
      q$right[current_q] <- mass[gap_idx]
      current_q <- current_q + 1L
    }
    #condition for the ending point of the right
    if (condition_num(mass = mass[(gap_idx + 1L):right], times, sample_num)) {
      current_boundary <- current_boundary + 1L
      boundary$left[current_boundary] <- gap_idx + 1L
      boundary$right[current_boundary] <- right
    }
    else {
      #mass[(gap_idx + 1L):right] <- r
      q$left[current_q] <- mass[gap_idx]
      q$right[current_q] <- mass[right]
      current_q <- current_q + 1L
    }
  }

  q$left <- q$left[1:current_q - 1L]
  q$right <- q$right[1:current_q - 1L]
  return(q)
}

#' Binning function
#'
#' @param mass
#' @param intensities
#' @param samples
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
binning_2 <- function(mass, intensities, samples, tolerance) {
  n <- length(mass)
  d <- diff(mass)
  n_boundaries <- max(20L, floor(3L * log(n)))
  boundary <- list(left = double(n_boundaries),
                   right = double(n_boundaries))
  current_boundary <- 1L
  boundary$left[current_boundary] <- 1L
  boundary$right[current_boundary] <- n
  while (current_boundary > 0L) {

    left <- boundary$left[current_boundary]
    right <- boundary$right[current_boundary]
    current_boundary <- current_boundary - 1L
    gaps <- d[left:(right - 1L)]
    gap_idx <- which.max(gaps) + left - 1L

    l <- condition(mass = mass[left:gap_idx],
                   samples = samples[left:gap_idx],
                   tolerance = tolerance)
    if (is.na(l)) {
      current_boundary <- current_boundary + 1L
      boundary$left[current_boundary] <- left
      boundary$right[current_boundary] <- gap_idx
    }
    else {
      mass[left:gap_idx] <- l
    }

    r <- condition(mass = mass[(gap_idx + 1L):right],
                   samples = samples[(gap_idx + 1L):right],
                   tolerance = tolerance)
    if (is.na(r[1L])) {
      current_boundary <- current_boundary + 1L
      boundary$left[current_boundary] <- gap_idx + 1L
      boundary$right[current_boundary] <- right
    }
    else {
      mass[(gap_idx + 1L):right] <- r
    }
    if (current_boundary == n_boundaries) {
      n_boundaries <- floor(n_boundaries * 1.5)
      boundary$left <- c(boundary$left,
                         double(n_boundaries - current_boundary))
      boundary$right <- c(boundary$right,
                          double(n_boundaries - current_boundary))
    }
  }
  return(mass)
}

#' Bin peaks
#'
#' @param l
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
binPeaks <- function(l, tolerance = 0.002) {
  nn <- sapply(l, nrow)
  non_empty <- nn != 0L
  samples <- rep.int(seq_along(l), nn)
  mass <- unname(unlist((lapply(l[non_empty], function(x) x$mz)),
                        recursive = FALSE, use.names = FALSE))
  intensities <- unlist(lapply(l[non_empty], function(x) x$intensity),
                        recursive = FALSE, use.names = FALSE)
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  mass <- binning_2(mass = mass, intensities = intensities, samples = samples,
                    tolerance = tolerance)
  s <- sort.int(mass, index.return = TRUE)
  mass <- s$x
  intensities <- intensities[s$ix]
  samples <- samples[s$ix]
  lIdx <- split(seq_along(mass), samples)
  l[non_empty] <- mapply(FUN = function(p, i) {
    p$mz <- mass[i]
    p$intensity <- intensities[i]
    p
  }, p = l[non_empty], i = lIdx, MoreArgs = NULL, SIMPLIFY = FALSE,
  USE.NAMES = FALSE)
  l
}

#' Binning condition
#'
#' @param mass
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
condition_ppm <- function(mass, tolerance) {
  ## If the ppm between the first and the last mz of the bin is larger than
  ## the threshold ---> continue binning
  l <- length(mass)
  if (abs((mass[l] - mass[1L]) / mass[l]) * 1e6 > tolerance) {
    return(TRUE)
  }
  FALSE
}

#' Title
#'
#' @description if the number of mz in this bin is larger than the threshold
#' ---> continue binning
#' the floor of the times is 1 (ideally the lowest peak number in a bin =
#' the sample number )
#' the cell of the times depends on the RAM of your computer and your
#' expectation of running speed
#'
#' @param mass
#' @param times
#' @param sample_num
#'
#' @return
#' @export
#'
#' @examples
condition_num <- function(mass, times, sample_num) {
  l <- length(mass)
  if (l > sample_num * times) {
    return(TRUE)
  }
  return(FALSE)
}
