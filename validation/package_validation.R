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

# *** 2.1 Extract and align scans -----------------------------------------
# extract the spectra information from the mzml files using default values for 
# binning, low intensity, and missingness (per measurement)
#qc_files <- list.files(".",pattern = "sumr_qc.*mzML")
#qc <- sumR::mzml_extract_magic(qc_files, cores = 6)

mzL <- sumR::mzml_extract_magic(cores = 6)
save(mzL, file = "mzL.Rdata")
load("mzL.Rdata")

# create positive and negative mz object data frames from the list of spectra
mz_obj_n <- mzL$neg
# 266,678 peaks
mz_obj_p <- mzL$pos
# 285,312

# *** 2.2 Technical quality check -----------------------------------------
tech_quality_n <- mz_quality_magic(mz_obj_n, meta_neg)
tech_quality_p <- mz_quality_magic(mz_obj_p, meta_pos)


# *** 2.3 Blank Subtraction -----------------------------------------------
mz_obj_n <- mz_subtraction(mz_obj_n, 
                           dplyr::filter(meta_neg, type != "IS_Blank" & 
                                           filtered_out == "no"),
                           dplyr::filter(meta_neg, type == "IS_Blank" &
                                           filtered_out == "no"))
# peaks before blank sub: 266,678
# peaks after blank sub: 45,218
# peaks removed by blank sub: 221,460

# *** 2.4 Mass defect -----------------------------------------------------
mz_obj_n <- mz_mass_defect(mz_obj_n, plot = TRUE)
# peaks before mass defect: 45,218
# peaks after mass defect: 15,220
# peaks removed by mass defect: 29,998

# *** 2.5 Missingness -----------------------------------------------------
mz_mdobj_n <- mz_filter_missingness(mz_obj_n, threshold = 0.1)
# missingness filter not validated due to issue with alignment/binning

# *** 2.6 Low intensity ---------------------------------------------------
mz_liobj_n <- mz_filter_lowIntensity(mz_mdobj_n, threshold = 1200)
# after blank sub: 45,218
# after low intensity: 44,453
# peaks removed by low intensity: 765

# *** 2.6 Filter Magic check ---------------------------------------------------
mz_magic <- mz_filter_magic(mz_obj_n, min_intensity = 1200, 
                            missingness_threshold = 0.1, blacklist = F)

testthat::expect_identical(mz_magic, mz_liobj_n)

# 3. Post processing -----------------------------------------------------
# *** 3.1 Imputation ------------------------------------------------------
mz_impobj_n <- mz_pp_imputation(mz_obj_n, low_noise = 12, high_noise = 1000)

# *** 3.2 Normalization ---------------------------------------------------
mz_norm_n <- mz_pp_normalization(mz_impobj_n, metadata = meta_neg, plot = TRUE)
# issue with plotting, so normalization not validated yet

# *** 3.3 Log transform ---------------------------------------------------
mz_long_n <- mz_pp_pivot_longer(mz_norm_n, plot = FALSE)
mz_long_n <- mzlong_pp_log_transform(mz_long_n, plot = FALSE)
mz_log_n <- mz_pp_log(mz_norm_n)

# *** 3.4 PostProcess Magic check ---------------------------------------------------
mz_magic <- mz_pp_magic(mz_impobj_n, meta_neg, noise = c(12,1000), plot = F)

testthat::expect_identical(mz_magic,list("mzLong" = mz_long_n ,"mzLog" = mz_log_n)
)

# 4. Statistical Analysis -------------------------------------------------

# *** 4.1 PCA -------------------------------------------------------------
mzlog_analysis_pca(mz_log_n, meta_neg)

# *** 4.2 T-Test ----------------------------------------------------------
welch_n <- mzlog_analysis_welch(mz_log_n,  
                     dplyr::filter(meta_neg, phenotype == "HepG2" & 
                                   filtered_out == "no"),
                     dplyr::filter(meta_neg, phenotype == "HEK293T" &
                                   filtered_out == "no"))

# *** 4.3 Fold & Volcano ----------------------------------------------------------
fold_n <- mzlog_analysis_fold(mz_log_n,  
                    dplyr::filter(meta_neg, phenotype == "HepG2" & 
                                    filtered_out == "no"),
                    dplyr::filter(meta_neg, phenotype == "HEK293T" &
                                    filtered_out == "no"))

plot_volcano(welch_n, fold_n, title = "Volcano Plot")
  
# *** 4.4 RandomForest ----------------------------------------------------------
