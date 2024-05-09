# ============
# identify_hmdb(mzs, adducts, hmdb_file)
# ============
test_that("identfy_hmdb: with known adducts", {
  mz <- c(295.1448, 254.4548, 60.89)
  ad <- c("M-H","3M+Na")
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz, ad)
	expect_identical(actual, hmdb_filtered)
})


test_that("identfy_hmdb: with custom adducts", {
  mz <- c(194.1440, 262.4190, 43.123)
  ad <- c(-100.0001, +0.301)
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz, ad)
	expect_identical(actual, hmdb_filtered)
})


# ============
# identify_lipids(mzs, adducts, hmdb_file)
# ============
test_that("identfy_lipids: with default adduct(H)", {
  mz <- c(174.0516, 497.4431, 199.199)
  ad <- c("HCOO-", "M+H")
  load("../test_data/asserts/lipid_filtered.Rdata")
  actual <- identify_lipids(mz,ad)
	expect_identical(actual, lipid_filtered)
})

# ============
# identify.adducts
# ============
test_that("identify.adducts: with named adducts", {
  ad <- c("HCOO-", "M+H")
  actual <- identify.adducts(ad,"testthat")
  expected <- c( -1.0008, 46.0254)
  expect_identical(actual, expected)
})

test_that("identify.adducts: with given mass", {
  ad <- c(123.123, 1701.1701)
  actual <- identify.adducts(ad,"testthat")
  expect_identical(actual, ad)
})

test_that("identify.adducts: with single value", {
  ad <- 100.012
  expect_identical(identify.adducts(ad,"testthat"), ad)
})

