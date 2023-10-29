# ============
# identify_hmdb(mzs, adducts, hmdb_file)
# ============
test_that("identfy_hmdb: with known adducts", {
  mz <- c(295.1449, 301.2163)
  ad <- c("H","K")
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz,ad)
	expect_identical(actual, hmdb_filtered)
})


test_that("identfy_hmdb: with custom adducts", {
  mz <- c(394.1440, 262.4181)
  ad <- c(100.0001, 0.3001)
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz,ad)
	expect_identical(actual, hmdb_filtered)
})
