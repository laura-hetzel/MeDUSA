###Reference zz_centroiding for the original/more information

#Return MZ,I list for single scan
centroid.singleScan <- function(spectrum, halfWindowSize = 2) {
  if (nrow(spectrum) == 0) return(NULL)

  intensities <- centroid.savgol(spectrum[, 2], halfWindowSize)
  intensities[intensities < 0] <- 0
  ## find local maxima
  centroids <- which(diff(sign(diff(intensities))) == -2) + 1

  spectrum <- as.data.frame(spectrum[centroids, ])
  colnames(spectrum) <- c("mz", "intensity")
  spectrum
}


centroid.savgol <- function(y, halfWindowSize = 2L, polynomialOrder = 3L) {

  sav.filter <- function(x, hws, coef) {
    n <- length(x)
    w <- 2L * hws + 1L
    y <- stats::filter(x = x, filter = coef[hws + 1L, ], sides = 2L)
    attributes(y) <- NULL
    y[seq_len(hws)] <- head(coef, hws) %*% head(x, w)
    y[(n - hws + 1L):n] <- tail(coef, hws) %*% tail(x, w)
    y
  }

  coef <- function(m, k = 3L) {
    k <- 0L:k
    nm <- 2L * m + 1L
    nk <- length(k)
    K <- matrix(k, nrow = nm, ncol = nk, byrow = TRUE)
    F <- matrix(double(), nrow = nm, ncol = nm)
    for (i in seq_len(m + 1L)) {
      M <- matrix(seq_len(nm) - i, nrow = nm, ncol = nk, byrow = FALSE)
      X <- M^K
      T <- solve(t(X) %*% X) %*% t(X)
      F[i, ] <- T[1L, ]
    }
    F[(m + 2L):nm, ] <- rev(F[seq_len(m), ])
    F
  }

  sav.filter(y, hws = halfWindowSize, coef = coef(m = halfWindowSize,
                                                  k = polynomialOrder))
}
