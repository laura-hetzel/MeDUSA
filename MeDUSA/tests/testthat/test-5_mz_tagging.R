##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_tag_isotope_hunter                                                        ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_tag_isotope_hunter: C13 Happy Path ", {
  load("testdata/mz_neg.Rdata")
  test_subset<-round(mz_neg[seq(240,680,10),],4)
  expect <- setNames(c(189.0771,193.0905,506.3499,507.3540),c(189.07711,189.07712,506.34991,506.34992))
  actual <- unlist(mz_tag_isotope_hunter(test_subset))
  expect_identical(actual,expect)
})

test_that("mz_tag_isotope_hunter: C13 Happy Path; only MZ", {
  load("testdata/mz_neg.Rdata")
  test_subset<-round(mz_neg[seq(240,680,10),],4)
  expect <- setNames(c(189.0771,193.0905,506.3499,507.3540),c(189.07711,189.07712,506.34991,506.34992))
  actual <- unlist(mz_tag_isotope_hunter(test_subset$mz))
  expect_identical(actual,expect)
})

test_that("mz_tag_isotope_hunter: custom params", {
  load("testdata/mz_neg.Rdata")
  test_subset<-round(mz_neg[seq(240,680,10),],4)
  expect <- setNames(c(205.1600, 214.1533,227.0387, 245.0244),c(205.161,205.162,227.03871,227.03872))
  actual <- unlist(mz_tag_isotope_hunter(test_subset, iso_target=2.9979, tol_ppm = 5e-5, iso_iter = 6))
  expect_identical(actual,expect)
})
