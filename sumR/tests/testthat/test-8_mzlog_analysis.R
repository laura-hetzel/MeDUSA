test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


# ============
# mzlog_analysis_pca(input_mzlog_obj,metadata, sample_blacklist)
# ============
test_that("mzlog_analysis_pca: HappyPath", {
  mz <- ##create/find test-data
  meta <- ##create/find test-data
  expected <- "wat" ##create/find test-data
  actual <- "wat"   #mzlog_analysis_pca(mz,meta)
  expect_identical(actual, expected)
})
