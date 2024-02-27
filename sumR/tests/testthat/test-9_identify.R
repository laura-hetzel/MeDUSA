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
