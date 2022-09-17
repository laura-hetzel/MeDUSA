library(sumR)
library(testthat)
files <- list.files(system.file("cells", package = "sumR"), full.names = TRUE)

test_that("mzML files can be centroided", {
  peaks <- extractPeaks(files)
  expect_equal(typeof(peaks), "list")
  expect_equal(length(peaks), length(files))

  expect_equal(class(peaks[[1]]), "data.frame")
  expect_true(all(c('scan', "rt", "mz", "i") %in% colnames(peaks[[1]])))
  attrs <- names(attributes(peaks[[1]]))
  expect_true(all(c("polarity", "massWindow", "combineSpectra", "Files", "rt",
    "scanranges", "Datetime") %in% attrs))
})

test_that("signals can be aligned across spectra", {

})

test_that("Peaks can be aligned across cells", {

})

test_that("Imputation can be done" , {
  # imputation(method = "noise", noise = 100, seed = 42)
})

test_that("Blank substraction is done correctly", {
   # blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE)
})

test_that("Mass Defect filter works properly", {
  # massDefectFilter()
})

test_that("Isotope tagging works", {
  # isotopeTagging(corr = 0.8)
})

test_that("Fragments are filtered", {
  # fragmentFilter(method = "spearman", corr = 0.95)
})
