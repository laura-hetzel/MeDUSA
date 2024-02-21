# ============
# mz_subtraction((sample_mz_obj, subtract_mz_obj , method = mean, threshold = 3)
# ============
testthat("mz_subtraction: Happy path",{
  load("../test_data/mz_neg.Rdata")
  load("../test_data/asserts/mz_subtraction.Rdata")
  meta <- data.frame(read.csv("../test_data/meta_neg.csv"))
  actual <- mz_subtraction(mztools_filter(mz_neg, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_neg, meta,"IS_Blank","type",F,T))
                           
  expect_identical(actual,mz_subtraction)
  
})

