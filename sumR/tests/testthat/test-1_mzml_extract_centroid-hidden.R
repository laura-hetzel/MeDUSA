# ============
# polarity_loop / polarity_final
# ============
# Highly coupled to extract_magic.

# ============
# fill_defaults(params)
# ============

test_that("fill_defaults: merges correctly",{
  expected <- list(
    "prebin_method" = acos,
    "intensity_threshold" = 1,
    "postbin_method" = max,
    "tolerance" = 5e-6,
    "timeSquash_method" = mean,
    "missingness_threshold" = .1,
    "massWindow" = c(0, Inf)
  )

  in_params <- list(
    "prebin_method" = acos,
    "intensity_threshold" = 1
  )
  actual <- magic.fill_defaults(in_params)
  expect_identical(actual,expected)
})

test_that("fill_defaults: fully replaces if needed",{
  expected <- list(
    "prebin_method" = max,
    "postbin_method" = max,
    "tolerance" = 5e-6,
    "timeSquash_method" = mean,
    "missingness_threshold" = .1,
    "intensity_threshold" = 1000,
    "massWindow" = c(0, Inf)
  )
  actual <- magic.fill_defaults()
  expect_identical(actual,expected)
})

# ============
# magic.binning(df, method, log_name, tolerance = 5e-6)
# ============
test_that("binning: handles custom params",{
  df <- data.frame( "mz"   = c(50, 51, 55, 200, 1000),
                    "sam1" = c(10, 20, 0,  0,   0),
                    "sam2" = c(0,  30, 40, 50,  0),
                    "sam3" = c(0,  0,  60, 70,  0),
                    "sam4" = c(0,  0,  0,  80,  90))
  expected <- data.frame( "mz"   = c(52,  200, 1000),
                          "sam1" = c(40,  0,   0),
                          "sam2" = c(80,  100, 0),
                          "sam3" = c(120, 140, 0),
                          "sam4" = c(0,   160, 180))
  max_double <- function(x, na.rm, na.action){ max(x) * 2 }

  actual <- magic.binning(df, max_double, tolerance = 0.5)
  expect_identical(actual,expected)
})

# ============
# magic.file_list(data)
# ============
test_that("file_list: lists from dir", {
  expect <- c("../test_data/mzml_1.mzML",  "../test_data/mzml_2.mzML",
              "../test_data/mzml_31.mzML", "../test_data/mzml_33.mzML")
  actual <- magic.file_lister("../test_data")
  expect_identical(actual,expect)
})

test_that("file_list: keeps file list", {
  expect <- list.files("../test_data", pattern = "mzml_3.*mzML", full.names=T)
  actual <- magic.file_lister(expect)
  expect_identical(actual,expect)
})

test_that("file_list: does not accept non-mzml files", {
  expect_error(
    magic.file_lister(c("mzml_1.llama", "mzml_2.chicken"),"ERROR: Cannot find mzML files in given files.")
  )
})

# ============
# magic.mz_format(data)
# ============
test_that("mz_format: flattens correctly", {
  input <- list(list(mz=c(1.1, 1.2, 5.3, 9.4), scan1=c(100.133252,200.1252352,300.5253252,400.853425)),
                list(mz=c(1.1, 3.5111111, 9.88888888), scan2=c(700.541345,800.321344,900.9243525)))
  expect <- data.frame( "mz" = c(1.1, 1.1, 1.2, 3.51111, 5.3, 9.4, 9.88889),
                     "scan1" = c(100.133252, NA, 200.1252352, NA ,300.5253252, 400.853425, NA),
                     "scan2" = c(NA, 700.541345 , NA ,800.321344, NA, NA, 900.9243525))
  expect_identical(magic.mz_format(input), expect)
})
