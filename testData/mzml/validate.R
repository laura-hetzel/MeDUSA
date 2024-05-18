library('MeDUSA')

# 2. Data Processing (filtering) ------------------------------------------
# set the working directory to where the mzml files are located
setwd("~/local/testData/mzml")
# introduce the meta data to use for filtering
file_meta <- "meta.xlsx"
meta <- data.frame(readxl::read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta$sample_name <- meta$filename
meta <- dplyr::select(meta,-filename)

###
# Extraction Options ( magic, MzT, load )
###
#Load existing
load("mzL.Rdata")
#Magic
#mzL <- mzml_extract_magic("data")
## debugging without magic
#files <- list.files(path=getwd("data"), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))

mz_obj_p <- mzL$pos

tech_quality_p <- mz_quality_magic(mz_obj_p, meta)

mz_obj_p <- mz_subtraction(mztools_filter(mz_obj_p, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_obj_p, meta,"IS_Blank","type",F,T))

mz_obj_p <- mz_mass_defect(mz_obj_p, plot = TRUE)
mz_obj_p <- mz_filter_magic(mz_obj_p, blacklist= F)
pp_out <- mz_post_magic(mz_obj_p, metadata = meta, plot = TRUE)

mzlog_analysis_pca(pp_out$mzLog, meta)

welch_p <- mzlog_analysis_welch( mztools_filter(pp_out$mzLog,meta,"HepG2"),
                                 mztools_filter(pp_out$mzLog,meta,"HEK293T"))

fold_p <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                              mztools_filter(pp_out$mzLog,meta,"HEK293T"))

plot_volcano(welch_p, fold_p)

rf_data <- mztools_filter(pp_out$mzLog, meta, c("HepG2","HEK293T"))

##RandomForest:Magic
rf_mag <- mzlog_rf_magic(rf_data, meta, "pos")
##RandomForest: without magic
#rf_cor <- mzlog_rf_correlation(rf_data)
#rf_sel <- mzlog_rf_select(rf_cor, meta)
#rf_obj <- mzlog_rf(rf_data, meta, rfe_obj = rf_sel)
#rf_validate <- rf_validate(rf_obj)
#rf_permuted <- rf_permuted(rf_obj, "pos")

anova <- mzlong_analysis_anova(pp_out$mzLong, meta, phenotypes=c("HepG2","HEK293T"))

hmdb_rf <- identify_hmdb(rf_mag$mz, c("M-H","M+Na")) 
hmdb_anova <- identify_hmdb(anova$imp_mz) 
lipids_rf <- identify_lipids(rf_mag$mz)
lipids_rf <- identify_lipids(anova$imp_mz)



