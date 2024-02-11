
# ============
# mztools_filter(input_mzobj, metadata, filter_value , filter_name = "phenotype", exclude = F, keep_mz = T)
# ============
test_that("mztools_filter: Filter By Phenotype", {
  load("../test_data/mz_neg.Rdata")
  load("../test_data/asserts/mztools_pheno_filter_red.Rdata")
  #meta_csv <- load("../test_data/meta_neg.csv")
  meta_csv <- read.csv("../test_data/meta_neg.csv")
  actual <- mztools_filter(mz_neg, meta_csv, "red")
  expect_identical(actual, phenotype_filter_red)
})

test_that("mztools_filter: Filter multiple Phenotypes", {
  load("../test_data/mz_neg.Rdata")
  load("../test_data/asserts/mztools_3phenos.Rdata")
  meta_csv3 <- read.csv("../test_data/meta_neg_3phenos.csv")
  actual <- mztools_filter(mz_neg, meta_csv3, c("red","blue"))
  expect_identical(actual, mztools_3phenos)
})


test_that("mztools_filter: Exclude By Phenotype", {
    load("../test_data/mz_neg.Rdata")
    load("../test_data/asserts/mztools_pheno_filter_red.Rdata")
    meta_csv <- read.csv("../test_data/meta_neg.csv")
    actual <- mztools_filter(mz_neg, meta_csv, "red", exclude = T)
    expect_identical(actual, phenotype_filter_red)
  })

test_that("mztools_filter: Drop MZ column", {
  load("../test_data/mz_neg.Rdata")
  meta_csv <- read.csv("../test_data/meta_neg.csv")
  load("../test_data/asserts/mztools_drop_mz.Rdata")
  actual <- mztools_filter(mz_neg, meta_csv, "red", keep_mz = F)
  expect_identical(actual, mztools_mz_col)
})

test_that("mztools_filter: Filter By Custom Column", {
  load("../test_data/mz_neg.Rdata")
  meta_csv <- read.csv("../test_data/meta_neg.csv")
  load("../test_data/asserts/mztools_custom_filter.Rdata")
  actual <- mztools_filter(mz_neg, meta_csv, 22, filter_name = "time")
  expect_identical(actual, mztools_custom_filter)
})

test_that("mztools_filter: [negative tests]", {
  load("../test_data/mz_neg.Rdata")
})

# use the below code as a reference
# test_that("mz_quality_metrics:  ", {
#   load("../test_data/mz_neg.Rdata")
#   load("../test_data/asserts/quality_check.Rdata")
#   actual <- mz_quality_metrics(mz_neg)
#   expect_identical(actual, quality_check)
# })
