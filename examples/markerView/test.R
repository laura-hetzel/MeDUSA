
library("MeDUSA")
setwd("~/local/examples/markerView")
data_file <-  "neg.csv"
meta_file <-  "meta.csv"
mz_obj <- data.frame(read.csv(data_file))
meta <- data.frame(read.csv(meta_file))
colnames(mz_obj) <- gsub("\\.txt$", "", colnames(mz_obj))

mz_obj <- tibble::column_to_rownames(mz_obj, "mz")
meta <- dplyr::rename(meta,c(sample_name=filename))

mz_obj <- mz_subtraction( mztools_filter(mz_obj, meta,"solvent","type",T,T),
                          mztools_filter(mz_obj, meta,"solvent","type",F,T))
mz_obj <- mz_subtraction( mztools_filter(mz_obj, meta,"cell","type",F,T),
                          mztools_filter(mz_obj, meta,"media_cell","type",F,T))
mz_obj <- mz_mass_defect(mz_obj, plot = TRUE)
mz_obj <- mz_filter_magic(mz_obj, blacklist= F)
pp_out <- mz_post_magic(mz_obj, metadata = meta, plot = TRUE)

bad_neg_samples <- c("sample155_seg_fullscan_neg",
                     "sample075_seg_fullscan_neg",
                     "sample168_seg_fullscan_neg",
                     "sample181_seg_fullscan_neg",
                     "sample164_seg_fullscan_neg")
mzlog_analysis_pca(mztools_filter(pp_out$mzLog,meta,"cell","type"), meta, sample_blacklist = bad_neg_samples)

mzlog_analysis_volcano_magic(mztools_filter(pp_out$mzLog,meta,"differentiated"),
                             mztools_filter(pp_out$mzLog,meta,"naive"))

welch <- mzlog_analysis_welch(mztools_filter(pp_out$mzLog,meta,"differentiated"),
                              mztools_filter(pp_out$mzLog,meta,"naive"))

fold <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"differentiated"),
                            mztools_filter(pp_out$mzLog,meta,"naive"))
plot_volcano(welch, fold)

rf_data <- mztools_filter(pp_out$mzLog, meta, c("differentiated","naive"))
rf_mag <- mzlog_rf_magic(rf_data, meta, "neg")
anova <- mzlong_analysis_anova(pp_out$mzLong, meta, phenotypes=c("differentiated","differentiated"))


hmdb_rf <- identify_hmdb(rf_mag$mz, c("M-H","M+Na"))
hmdb_anova <- identify_hmdb(anova$imp_mz)
lipids_rf <- identify_lipids(rf_mag$mz)
lipids_rf <- identify_lipids(anova$imp_mz)

