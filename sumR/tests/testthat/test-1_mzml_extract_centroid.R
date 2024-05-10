##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mzml_extract_magic                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzml_extract_magic: threaded, in current directory (all defaults)",{
  setwd("testdata")
  load("asserts/mzL_all.Rdata")
  actual <- mzml_extract_magic()
  expect_identical(actual,mzL_all)
})

test_that("mzml_extract_magic: single thread with file filter",{
  load("testdata/mzL_3.Rdata")
  files <- list.files("testdata", pattern = "mzml_3.*mzML")
  actual <- mzml_extract_magic(files,cores=1)
  expect_identical(actual,mzL_3)
})

test_that("mzml_extract_magic: handles custom param in another dir",{
  load("testdata/asserts/mzL_custom.Rdata")
  load("testdata/asserts/mzL_all.Rdata")
  actual <- mzml_extract_magic("testdata",cores=4,params=c("intensity_threshold" = 100000 ))
  expect_identical(actual,mzL_custom)
  expect_false(isTRUE(all.equal(mzL_all, mzL_custom)))
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mzml_extract_file                                                    ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mzml_extract_file: can extract all (no polarity filter) with no magic",{
  load("testdata/asserts/mzT_1.Rdata")
  actual <- mzml_extract_filec("testdata/mzml_1.mzML")
  expect_identical(actual,mzT_1)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mzT_squashTime                                                    ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO: test_that("mzT_squashTime: with custom method"{})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mzT_filtering                                                    ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Adequately tested by the above
