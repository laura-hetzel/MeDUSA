##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.mz_log_removed_rows                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_log_removed_rows: HappyPath", {
  inp <- data.frame( mz   <- c(50,100,150,200,1000),
                     sam1 <- c(10, 20, 0 , 0,  0),
                     sam2 <- c(0,  30, 40, 0,  0),
                     sam3 <- c(0,  0,  50, 60, 0),
                     sam4 <- c(0,  0,  80, 90, 100)
                   )
  expect_no_error( local.mz_log_removed_rows(inp,inp[1,3],"Clever Accounting") )
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.mz_polarity_guesser                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_polarity_guesser: HappyPath",  {
    load("testdata/mz_neg.Rdata")
    expect_identical(local.mz_polarity_guesser(mz_neg), "Negative")
    colnames(mz_neg) <- c('mz',"Positive_sample1","Positive_sample2","Positive_sample3", "Positive_sample4")
    expect_identical(local.mz_polarity_guesser(mz_neg), "Positive")
})

test_that("mz_polarity_guesser: Custom Return",  {
    load("testdata/mz_neg.Rdata")
    expect_identical(local.mz_polarity_guesser(mz_neg, neg_return = "Chickens"), "Chickens")
})

test_that("mz_polarity_guesser: Errors on both",  {
    load("testdata/mz_neg.Rdata")
    colnames(mz_neg) <- c('mz',"sample1","sample2","sample3", "sample4")
    expect_error(local.mz_polarity_guesser(mz_neg), "ERROR: MeDUSA::polarity_guesser: Could not guess Positive or Negative from colnames")
})

test_that("mz_polarity_guesser: Errors on both",  {
    load("testdata/mz_neg.Rdata")
    colnames(mz_neg) <- c('mz',"neg_sample1","Positive_sample2","Positive_sample3", "Positive_sample4")
    expect_error(local.mz_polarity_guesser(mz_neg), "ERROR: MeDUSA::polarity_guesser: Detected both Positive & Negative from colnames")
})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.meta_polarity_fixer                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_polarity_guesser: HappyPath",  {
    load("testdata/mz_neg.Rdata")
    meta <- data.frame(read.csv("testdata/meta_neg.csv"))
    expect <- c("neg_mzml_1", "neg_mzml_2",  "neg_mzml_31", "neg_mzml_33")
    expect_identical(local.meta_polarity_fixer(mz_neg, meta)$sample_name, expect)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.ensure_mz                                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("ensure_mz: no_mz provided",{
    load("testdata/mz_neg.Rdata")
    meta <- data.frame(read.csv("testdata/meta_neg.csv"))
    expect_error(
        local.ensure_mz(
            mztools_filter(mz_neg, meta, c("blue"), keep_mz = F),
            mztools_filter(mz_neg, meta, c("red"),  keep_mz = F), "llama" ),
        "ERROR: llama: neither input dataframes have column mz")
})

test_that("ensure_mz: mismatch rows",{
    load("testdata/mz_neg.Rdata")
    meta <- data.frame(read.csv("testdata/meta_neg.csv"))
    expect_error(
        local.ensure_mz(
            mztools_filter(mz_neg[1:10,], meta, c("blue")),
            mztools_filter(mz_neg[1:9,], meta, c("red")), "llama" ),
        "ERROR: llama: Dataframes have a different number of mz")
})

test_that("ensure_mz: happyPath",{
    load("testdata/mz_neg.Rdata")
    load("testdata/asserts/ensure_mz_happy.Rdata")
    meta <- data.frame(read.csv("testdata/meta_neg.csv"))
    actual <-  local.ensure_mz(
        mztools_filter(mz_neg, meta, c("blue")),
        mztools_filter(mz_neg, meta, c("red")))
    expect_identical( actual, ensure_mz_happy )
})
