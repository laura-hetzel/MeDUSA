# ============
# identify_hmdb(mzs, adducts, hmdb_file)
# ============
test_that("identfy_hmdb: positive with known adducts", {
  mz <- c(295.1449, 301.2163)
  ad <- c("H","K")
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz, ad, "neg")
	expect_identical(actual, hmdb_filtered)
})


test_that("identfy_hmdb: negative with custom adducts", {
  mz <- c(194.1440, 261.8181)
  ad <- c(100.0001, 0.3001)
  load("../test_data/asserts/hmdb_filtered.Rdata")
  actual <- identify_hmdb(mz, ad, "pos" )
	expect_identical(actual, hmdb_filtered)
})


# ============
# identify_lipids(mzs, adducts, hmdb_file)
# ============
test_that("identfy_lipids: with default adduct(H)", {
  mz <- c(872.5805,501.2620)
  load("../test_data/asserts/lipid_filtered.Rdata")
  actual <-identify_lipids(mz)
  actual <- identify_hmdb(mz,ad)
	expect_identical(actual, lipid_filtered)
})
