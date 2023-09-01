library("sqrlSumr")

mzL <- sqrlSumr::mzml_extract_magic(cores=4)

meta <- data.frame(read.csv("meta_neg.csv"))

mz_obj <-  mzL$neg

#mz_carbon <- sqrlSumr::mz_isotope_hunter(mz_obj, cores=4)

mz_obj <- sqrlSumr::mz_subtraction(mz_obj,
            dplyr::filter(meta, type == "cell" & filtered_out == "no"),
            dplyr::filter(meta, type == "blank" & filtered_out == "no"))

mz_obj <- sqrlSumr::mz_mass_defect(mz_obj)

#misssingness, LowIntensity
mz_obj <- sqrlSumr::mz_filter_magic(mz_obj, blacklist=c(100,110))

#post processing
mz_mag <- sqrlSumr::mz_post_process_magic(mz_obj,meta)
mz_mag <- sqrlSumr::mz_pp_imputation(mz_obj)
mz_mag <- sqrlSumr::mz_pp_normalization(mz_mag, meta)
mz_mag <- log2(mz_mag)


## debugging without magic
#files <- list.files(path=getwd(), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
