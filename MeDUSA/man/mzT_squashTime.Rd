% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_mzml_extract_centroid.R
\name{mzT_squashTime}
\alias{mzT_squashTime}
\title{Squash MzT obj (many scans) into a single mz|intensity dataframe}
\usage{
mzT_squashTime(mzT, timeSquash_method = mean, ignore_zeros = T)
}
\arguments{
\item{mzT}{\cr
MzT : MzT object}

\item{ignore_zeros}{\cr
Boolean: Should we set 0 <- NA (to not affect the math)}

\item{method}{\cr
\link{Math} : i.e. (mean, max, median)}
}
\description{
Squash MzT obj (many scans) into a single mz|intensity dataframe
}
