##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.mz_log_removed_rows                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_log_removed_rows: HappyPath", {
  in <-data.frame( mz   <- c(50,100,150,200,1000),
                   sam1 <- c(10, 20, 0 , 0,  0),
                   sam2 <- c(0,  30, 40, 0,  0),
                   sam3 <- c(0,  0,  50, 60, 0),
                   sam4 <- c(0,  0,  80, 90, 100),
                   )
  expect_no_error( local.mz_log_removed_rows(in,in[1,3],"Clever Accounting")
})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.mz_polarity_guesser                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_polarity_guesser", {

})

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.meta_polarity_fixer                                                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
local.meta_polarity_fixer <- function(input_mz, mesta){
  prepend <- local.mz_polarity_guesser(input_mz, pos_return = "pos", neg_return = "neg")
  meta$sample_name <- paste(prepend, meta$sample_name, sep="_")
  meta
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     local.ensure_mz                                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("ensure_mz: happy path",{


})

test_that("ensure_mz: mismatch rows",{

})

test_that("enxure_mz: no_mz pfovided",{

})
