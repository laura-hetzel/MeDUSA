library('sumR')

# 2. Data Processing (filtering) ------------------------------------------
# set the working directory to where the mzml files are located
setwd("~/local/testData/mzml")
# introduce the meta data to use for filtering
file_meta <- "meta.xlsx"
meta <- data.frame(readxl::read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta$sample_name <- meta$filename
meta <- dplyr::select(meta,-filename)

load("mzL.Rdata")

mz_obj_p <- mzL$pos

tech_quality_p <- mz_quality_magic(mz_obj_p, meta)

mz_obj_p <- mz_subtraction(mztools_filter(mz_obj_p, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_obj_p, meta,"IS_Blank","type",F,T))

mz_obj_p <- mz_mass_defect(mz_obj_p, plot = TRUE)
mz_obj_p <- mz_filter_magic(mz_obj_p, blacklist= F)
pp_out <- mz_pp_magic(mz_obj_p, metadata = meta, plot = TRUE)

mzlog_analysis_pca(pp_out$mzLog, meta)

welch_p <- mzlog_analysis_welch( mztools_filter(pp_out$mzLog,meta,"HepG2"),
                                 mztools_filter(pp_out$mzLog,meta,"HEK293T"))

fold_p <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                              mztools_filter(pp_out$mzLog,meta,"HEK293T"))

plot_volcano(welch_p, fold_p)

rf_data <- mztools_filter(pp_out$mzLog, meta, c("HepG2","HEK293T"))
rf_cor <- mzlog_rf_correlation(rf_data)
rf_sel <- mzlog_rf_select(rf_cor, meta)
rf_obj <- mzlog_rf(rf_data, meta, rfe_obj = rf_sel)
rf_validate <- rf_validate(rf_obj)
rf_permuted <- rf_permuted(rf_obj, "pos")

rf_mag <- mzlog_rf_magic(rf_data, meta, "pos")

anova <- mzlong_analysis_anova(pp_out$mzLong, meta, phenotypes=c("HepG2","HEK293T"))

