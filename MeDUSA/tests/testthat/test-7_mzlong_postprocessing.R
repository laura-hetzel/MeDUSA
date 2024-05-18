##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mzlong_pp_log_transform                                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzlong_pp_log_transform: HappyPath", {
  load("testdata/asserts/long_log_transform_happy.Rdata")
  load("testdata/asserts/pivot_longer_happy.Rdata")
  expect_identical(mzlong_pp_log_transform(pivot_longer_happy, F), long_log_transform_happy)
})
