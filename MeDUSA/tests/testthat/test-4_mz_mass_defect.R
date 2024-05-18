##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##     mz_mass_defect                                                       ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_that("mz_mass_defect: Happy Path", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mass_defect_happy.Rdata")
  actual <- mz_mass_defect(mz_neg,F)
  expect_identical(actual,mass_defect_happy)
})

test_that("mz_mass_defect: custom params", {
  load("testdata/mz_neg.Rdata")
  load("testdata/asserts/mass_defect_custom.Rdata")
  actual <- mz_mass_defect(mz_neg,F, 0.001,0.1)
  expect_identical(actual,mass_defect_custom)
})
