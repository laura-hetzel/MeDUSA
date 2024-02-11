# ============
# mz_filter_blacklist
# ============
test_that("mz_filter_blacklist: HappyPath", {
  load("../test_data/mz_neg.Rdata")
  load("../test_data/asserts/mz_filter_blacklist_happy.Rdata")
  actual <- mz_filter_blacklist(mz_neg[0:10,], c(56.99589,59.01289,59.02517, 1000))
  expect_identical(mz_filter_blacklist_happy, actual)
})

test_that("mz_filter_blacklist: Disallow Negative input", {
  load("../test_data/mz_neg.Rdata")
  expect_error( mz_filter_blacklist(mz_neg[0:10,],c(-90,1000000)) )
})

test_that("mz_filter_blacklist: happy blacklist file", {
  load("../test_data/mz_neg.Rdata")
  actual <- mz_filter_blacklist(mz_neg[0:10,],"../test_data/filter_blacklist_happy.csv")
  expect_identical(mz_filter_blacklist_happy, actual)
})

test_that("mz_filter_blacklist: bad data from file", {
  load("../test_data/mz_neg.Rdata")
  expect_error(mz_filter_blacklist(mz_neg[0:10,],"../test_data/filter_blacklist_badData.csv"))
})

test_that("mz_filter_blacklist: bad title from file", {
  load("../test_data/mz_neg.Rdata")
  expect_error(mz_filter_blacklist(mz_neg[0:10,],"../test_data/filter_blacklist_badTitle.csv"))
})

# ============
# mz_low_intensity
# ============

# ============
# mz_missingness
# ============

# ============
# mz_filter_magic
# ============
