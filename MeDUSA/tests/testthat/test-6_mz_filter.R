##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_filter_blacklist                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mz_filter_blacklist: HappyPath", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_blacklist_happy.Rdata")
  actual <- mz_filter_blacklist(mz_neg[0:10,], c(56.99589, 59.01289, 59.02517, 1000))
  expect_identical(mz_filter_blacklist_happy, actual)
})

test_that("mz_filter_blacklist: Disallow Negative input", {
  load("testdata/mz_neg.Rdata")
  expect_error( mz_filter_blacklist(mz_neg[0:10,], c(-90,1000000)) )
})

test_that("mz_filter_blacklist: happy blacklist file", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_blacklist_happy.Rdata")
  actual <- mz_filter_blacklist(mz_neg[0:10,], "testdata/filter_blacklist_happy.csv")
  expect_identical(mz_filter_blacklist_happy, actual)
})

test_that("mz_filter_blacklist: bad data from file", {
  load("testdata/mz_neg.Rdata")
  expect_error(mz_filter_blacklist(mz_neg[0:10,],"testdata/filter_blacklist_badData.csv"))
})

test_that("mz_filter_blacklist: bad title from file", {
  load("testdata/mz_neg.Rdata")
  expect_error(mz_filter_blacklist(mz_neg[0:10,],"testdata/filter_blacklist_badTitle.csv"))
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_filter_lowIntensity                                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mz_low_intensity: happy path", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_lowIntensity.Rdata")
  actual <- mz_filter_lowIntensity(mz_neg, 2500)
  expect_identical(actual,mz_filter_lowIntensity)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_filter_missingness                                                ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_missingness: percentage", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_missingness.Rdata")
  actual <- mz_filter_missingness(mz_neg, 0.5)
  expect_identical(actual,mz_filter_missingness)
})

test_that("mz_missingness: count", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_missingness.Rdata")
  actual <- mz_filter_missingness(mz_neg, 2)
  expect_identical(actual,mz_filter_missingness)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_filter_magic                                                      ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_filter_magic: defaults", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_magic.Rdata")
  actual <- mz_filter_magic(mz_neg)
  expect_identical(actual,mz_filter_magic)
})

test_that("mz_filter_magic: parameters", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_filter_magic_param.Rdata")
  actual <- mz_filter_magic(mz_neg, missingness_threshold=2, min_intensity=10000,
                            blacklist = c(56.99589,59.01289,59.02517, 1000))
  expect_identical(actual,mz_filter_magic_param)
})
