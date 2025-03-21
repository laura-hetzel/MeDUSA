% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_mzml_extract_centroid.R
\name{mzml_extract_magic}
\alias{mzml_extract_magic}
\title{Create MzObj from directory or file
Current Defaults:
mzBinningMath: max
mzBinTolerance = 5e-6
timeSquashMath: mean
missingness_threshold: 0.02 (2\%)
massWindow: c(0,Inf)}
\usage{
mzml_extract_magic(files = getwd(), cores = 2, ...)
}
\arguments{
\item{files}{\cr
String: File or directory of mzmL files}

\item{cores}{Integer: Can I has multithreading? (Need parallel)}
}
\description{
Create MzObj from directory or file
Current Defaults:
mzBinningMath: max
mzBinTolerance = 5e-6
timeSquashMath: mean
missingness_threshold: 0.02 (2\%)
massWindow: c(0,Inf)
}
