setwd("~/local/testData/mzml/1")
library("sumR")

mzL <- sumR::mzml_extract_magic(cores=4)

meta <- data.frame(read.csv("meta_neg.csv"))

mz_obj <-  mzL$neg

sumR::mz_quality_magic(mz_obj, meta, cores = 2)

#mz_carbon <- sumR::mz_isotope_hunter(mz_obj, cores=4)

mz_obj <- sumR::mz_subtraction(mz_obj,
            dplyr::filter(meta, type == "cell" & filtered_out == "no"),
            dplyr::filter(meta, type == "blank" & filtered_out == "no"))

mz_obj <- sumR::mz_mass_defect(mz_obj)

#misssingness, LowIntensity
mz_obj <- sumR::mz_filter_magic(mz_obj, blacklist=c(100,110))

#post processing
mz_mag <- sumR::mz_post_process_magic(mz_obj,meta)x
#mz_mag <- sumR::mz_pp_imputation(mz_obj)
#mz_mag <- sumR::mz_pp_normalization(mz_mag, meta)
#mz_mag <- log2(mz_mag)


## debugging without magic
#files <- list.files(path=getwd(), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
