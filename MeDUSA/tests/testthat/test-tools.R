##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mztools_filter                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mztools_filter: Filter By Phenotype", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mztools_pheno_filter_red.Rdata")
  #meta_csv <- load("testdata/meta_neg.csv")
  meta_csv <- read.csv("testdata/meta_neg.csv")
  actual <- mztools_filter(mz_neg, meta_csv, "red")
  expect_identical(actual, phenotype_filter_red)
})

test_that("mztools_filter: Filter multiple Phenotypes", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mztools_3phenos.Rdata")
  meta_csv3 <- read.csv("testdata/meta_neg_3phenos.csv")
  actual <- mztools_filter(mz_neg, meta_csv3, c("red","blue"))
  expect_identical(actual, mztools_3phenos)
})

test_that("mztools_filter: Exclude By Phenotype", {
    load("testdata/mz_neg.Rdata")
    load("testdata/asserts/mztools_pheno_filter_red.Rdata")
    meta_csv <- read.csv("testdata/meta_neg.csv")
    actual <- mztools_filter(mz_neg, meta_csv, "blue", exclude = T)
    expect_identical(actual, phenotype_filter_red)
  })

test_that("mztools_filter: Drop MZ column", {
  load("testdata/mz_neg.Rdata")
  meta_csv <- read.csv("testdata/meta_neg.csv")
  load("testdata/asserts/mztools_drop_mz.Rdata")
  actual <- mztools_filter(mz_neg, meta_csv, "red", keep_mz = F)
  expect_identical(actual, mztools_mz_col)
})

test_that("mztools_filter: Filter By Custom Column", {
  load("testdata/mz_neg.Rdata")
  meta_csv <- read.csv("testdata/meta_neg.csv")
  load("testdata/asserts/mztools_custom_filter.Rdata")
  actual <- mztools_filter(mz_neg, meta_csv, 22, filter_name = "time")
  expect_identical(actual, mztools_custom_filter)
})

#test_that("mztools_filter: [negative tests]", {
#  load("testdata/mz_neg.Rdata")
#})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     get_default_data                                                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("get_default_data: Mz Blacklist", {
  actual <- as.character(get_default_data('blacklist')[1,])
  expect_identical(actual, c("Polysiloxane", 536.17))
})

test_that("get_default_data: Adducts", {
  actual <- as.character(get_default_data('adducts')[3,])
  expect_identical(actual, c("M+NH4", -14.0067 ))
})
