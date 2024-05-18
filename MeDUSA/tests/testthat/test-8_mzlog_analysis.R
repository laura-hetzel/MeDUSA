##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      mzlog_analysis_pca                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mzlog_analysis_pca: HappyPath", {
  mz <- ##create/find test-data
  meta <- ##create/find test-data
  expected <- "wat" ##create/find test-data
  actual <- "wat"   #mzlog_analysis_pca(mz,meta)
  expect_identical(actual, expected)
})
