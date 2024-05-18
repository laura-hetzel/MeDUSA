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
test_that("mz_pp_normalization: HappyPath", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/normalization_happy.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  expect_identical(mz_pp_normalization(mz_neg, meta, F), normalization_happy)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_pivot_longer                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_pp_normalization: HappyPath", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/pivot_longer_happy.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  expect_identical(mz_pp_pivot_longer(mz_neg, F), pivot_longer_happy)
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_pp_magic                                                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_pp_normalization: HappyPath", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/pp_magic_happy.Rdata")
  meta <- data.frame(read.csv("testdata/meta_neg.csv"))
  expect_identical(mz_pp_magic(mz_neg, meta, noise = c(0,1), F), pp_magic_happy)
})
