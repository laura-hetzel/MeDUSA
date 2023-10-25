
# ============
# mzml_extract_magic(files, cores,  params )
#   downstream: mzml_extract_file, mzT_filtering, mzT_squashTime, magic.file_list
#               magic.binning, magic.mz_format, magic.polarity_loop,
#               magic.polarity_final, local.export_thread_env,
#               local.kill_threads, binning, centroid.singleScan,
#               mz_filter_lowIntensity, mz_filter_missingness
# ============

test_that("mzml_extract_magic: threaded, in current directory (all defaults)",{
  setwd("../test_data")
  load("asserts/mzL_all.Rdata")
  actual <- mzml_extract_magic()
  expect_identical(actual,mzL_all)
})

test_that("mzml_extract_magic: single thread with file filter",{
  load("../test_data/mzL_3.Rdata")
  files <- list.files("../test_data", pattern = "mzml_3.*mzML")
  actual <- mzml_extract_magic(files,cores=1)
  expect_identical(actual,mzL_3)
})

test_that("mzml_extract_magic: handles custom param in another dir",{
  load("../test_data/asserts/mzL_custom.Rdata")
  load("../test_data/asserts/mzL_all.Rdata")
  actual <- mzml_extract_magic("../test_data",cores=4,params=c("intensity_threshold" = 100000 ))
  expect_identical(actual,mzL_custom)
  expect_false(isTRUE(all.equal(mzL_all, mzL_custom)))
})

# ============
# mzml_extract_file(file, polarity, magic, cl,  params )
#   downstream: mzT_filtering, mzT_squashTime, magic.binning,
#               magic.mz_format,local.export_thread_env,
#               local.kill_threads, binning, centroid.singleScan
#               mz_filter_lowIntensity, mz_filter_missingness
# ============
test_that("mzml_extract_file: can extract all (no polarity filter) with no magic",{
  load("../test_data/asserts/mzT_1.Rdata")
  actual <- mzml_extract_filec("../test_data/mzml_1.mzML")
  expect_identical(actual,mzT_1)
})

# ============
# mzT_squashTime(mzT, timeSquash_method, ignore_zeros, cl)
# ============
# TODO: test_that("mzT_squashTime: with custom method"{})

# ============
# mzT_filtering(mzT, prebin_method, missingness_threshold, intensity_threshold, log_name)
#   downstream: magic.binning, mz_filter_lowIntensity, mz_filter_missingness
# ============
# Adequately tested by the above



