library("sqrlSumr")

meta <- data.frame(read.csv("../sum-r/testData/mzml/1/meta_neg.csv"))

mz_list <- mzml_extract_magic(cores=4)
mz_obj <- mz_list$neg

mz_obj <- sqrlSumr::mz_subtraction(mz_obj,
                                   dplyr::filter(meta, type == "blank" & filtered_out == "no"),
                                   dplyr::filter(meta, type == "cell" & filtered_out == "no"))

mz_obj <- sqrlSumr::mz_mass_defect(mz_obj)
mz_obj <- sqrlSumr::mz_filter_magic(mz_obj, blacklist=c(100,110))
