##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##      mzlong_analysis_anova                                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

test_that("mzlong_analysis_anova: HappyPath", {
    load("testdata/mz_neg.Rdata")
    meta <- data.frame(read.csv("testdata/meta_neg.csv"))
    expect <- c( 112.99)
    
    #TODO: make this independent of post processing
    mz_neg <- mz_post_imputation(mz_neg)
    mz_neg <- mz_post_pivot_longer(mz_neg)
    mz_neg <- mzlong_post_log(mz_neg)
    #TODO: this doesn't always produce the same result :/
    set.seed(105.9)
    actual <- mzlong_analysis_anova(mz_neg,meta,c('red','blue'))
    expect_identical(round(actual$imp_mz[1],2),expect)
})
