##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_quality_metrics                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_quality_metrics:  ", {
  load("testdata/mz_neg.Rdata")
	load("testdata/asserts/quality_check.Rdata")
  actual <- mz_quality_metrics(mz_neg)
	expect_identical(actual, quality_check)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_quality_meta_check                                                ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mz_quality_meta_check: HappyPath", {
  load("testdata/mz_neg.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  expect_no_error(mz_quality_meta_check(mz_neg,meta))
})

test_that("mz_quality_meta_check: Bad meta", {
  load("testdata/mz_neg.Rdata")
  meta <- data.frame(read.csv("testdata/asserts/meta_badneg.csv"))
  expect_error(mz_quality_meta_check(mz_neg,meta))
})
