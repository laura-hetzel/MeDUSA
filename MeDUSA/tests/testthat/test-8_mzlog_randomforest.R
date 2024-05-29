##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      mzlog_rf_correlation                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzlog_rf_correlation: Happy path, with param",{
  load("testdata/mz_neg.Rdata")
  expect <- c("sample_name", "123.0201775", "602.8247075")
  actual <- mzlog_rf_correlation(mz_neg,0.89)
  expect_identical(colnames(actual), expect )
})

test_that("mzlog_rf_correlation: Errors with too few correlated mz",{
  load("testdata/mz_neg.Rdata")
  expect_error(mzlog_rf_correlation(mz_neg))
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      mzlog_rf_select                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzlog_rf_select: Happy path",{
  load("testdata/mz_neg.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  #Variables Accuracy Kappa AccuracySD KappaSD
  expect <- c(2,0.3,0,0.47016,0)
  rf_cor <- mzlog_rf_correlation(mz_neg,0.89)
  actual <- round(mzlog_rf_select(rf_cor,meta)$results,5)
  expect_equal(sum( actual == expect), 5)
})


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      mzlog_rf                                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzlog_rf: Where are the trees?",{
  load("testdata/mz_neg.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  expect_error(mzlog_rf(mz_neg, meta_neg))
})
