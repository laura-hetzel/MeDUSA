##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_imputation                                                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Due to the random nature of imputation, it's asserts are difficult
test_that("mz_pp_imputation: Happy path", {
  load("testdata/mz_neg.Rdata")
  suppressWarnings(
    actual <- mz_pp_imputation(mz_neg[0:10,]) 
  )
  expect_lt( sum(actual<10) , 1 )
})

test_that("mz_pp_imputation: Parameters", {
  load("testdata/mz_neg.Rdata")
  actual <- mz_pp_imputation(mz_neg[0:10,], low_noise=150000, high_noise= 150001)
  expect_lt(sum(actual[colnames(actual) != 'mz'] < 130000 ), 1 )
  #Ensure MZ isn't imputated
  expect_identical( select(mz_neg[0:10,], mz) , select(actual,mz),)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_normalization                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_pivot_longer                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_magic                                                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
