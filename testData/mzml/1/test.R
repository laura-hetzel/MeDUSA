setwd("~/local/testData/mzml/1")
library("sumR")


# 1. General Info ---------------------------------------------------------
# The purpose of this script is to run MASDR (formerly known as SUMR) and compare
# the results to those obtained through scripts written for data processed in 
# MarkerView. This will validate MASDR as a processing tool in line with industry
# standards.


# *** 1.1 Sample Info -----------------------------------------------------
# cells - HEPG2 or Hek293T
# blank_solvent - sample of only the ionization solvent to serve as a baseline
#                 for which metabolites are in the solvent itself
# blank_media - media sampled from the either set of cells
# QC - quality control sample, diluted calibration mix

# *** 1.2 Dependencies ----------------------------------------------------
library('sumR')
library('readxl')
# 2. Data Processing (filtering) ------------------------------------------
# set the working directory to where the mzml files are located
setwd("~/local/mzml")
# introduce the meta data to use for filtering
file_meta <- "~/local/meta.xlsx"
meta_pos <- data.frame(read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta_pos$sample_name <- paste("pos",meta_pos$filename,sep="_")
meta_pos <- select(meta_pos,-filename)
meta_neg <- meta_pos
meta_neg$sample_name <- gsub("pos","neg",meta_neg$sample_name)
meta_neg$polarity <- "negative"
meta_pos$polarity <- "positive"


meta <- data.frame(read.csv("meta_neg.csv"))

mzL <- sumR::mzml_extract_magic(cores=4)
mz_obj <-  mzL$neg

sumR::mz_quality_magic(mz_obj, meta, cores = 2)


#mz_carbon <- sumR::mz_isotope_hunter(mz_obj, cores=4)

mz_obj <- sumR::mz_subtraction(mz_obj,
            dplyr::filter(meta, type == "cell" & filtered_out == "no"),
            dplyr::filter(meta, type == "blank" & filtered_out == "no"))

mz_obj <- sumR::mz_mass_defect(mz_obj)

#misssingness, LowIntensity
mz_obj <- sumR::mz_filter_magic(mz_obj, missingness_threshold = 3, blacklist=c(100,110))
#post processing
mz_mag <- sumR::mz_pp_magic(mz_obj,meta)
#RandomForest
cor_data <- mzlog_rf_correlation(mz_mag$mzLog, 0.9999)
rfe <- mzlog_rf_select(cor_data,meta,"type")

## debugging without magic
#files <- list.files(path=getwd(), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
