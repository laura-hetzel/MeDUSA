# ============
# mz_quality_metrics(input_mz_obj, cores)
# ============
test_that("mz_quality_metrics:  ", {
  load("../test_data/mz_neg.Rdata")
	load("../test_data/asserts/quality_check.Rdata")
  actual <- mz_quality_metrics(mz_neg)
	expect_identical(actual, quality_check)
})

# ============
# mz_quality_meta_check(input_mz_obj, meta)
# ============
test_that("mz_quality_meta_check: HappyPath", {
  load("../test_data/mz_neg.Rdata")
  meta <- data.frame(read.csv("../test_data/meta_neg.csv"))
  expect_no_error(mz_quality_meta_check(mz_neg,meta))
})

test_that("mz_quality_meta_check: Bad meta", {
  load("../test_data/mz_neg.Rdata")
  meta <- data.frame(read.csv("../test_data/asserts/meta_badneg.csv"))
  expect_error(mz_quality_meta_check(mz_neg,meta))
})
