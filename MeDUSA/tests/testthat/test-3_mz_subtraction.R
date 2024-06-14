# ============
# mz_subtraction((sample_mz_obj, subtract_mz_obj , method = mean, threshold = 3)
# ============
test_that("mz_subtraction: Happy path",{
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_subtraction.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  actual <- mz_subtraction(mztools_filter(mz_neg, meta,"red","phenotype",T,T),
                           mztools_filter(mz_neg, meta,"red","phenotype",F,T))
  expect_identical(actual,mz_subtraction)
})

test_that("mz_subtraction: Custom params",{
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mz_subtract_custom.Rdata")
  math <- function(x){mean(x)^22}
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  actual <- mz_subtraction(mztools_filter(mz_neg, meta,"red","phenotype",T,T),
                           mztools_filter(mz_neg, meta,"red","phenotype",F,T),
                           method = math, threshold = 2)

  expect_identical(actual,mz_subtract_custom)
})
