library('MeDUSA')

# Assumes Docker environment
# set the working directory to where the mzml files are located
# Don't forget to `git lfs pull` to get the example data
setwd("~/local/examples/mzml")
# introduce the meta data to use for filtering
file_meta <- "meta.xlsx"
meta <- data.frame(readxl::read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta$sample_name <- meta$filename
meta <- meta[meta$sample_name %in% 
       c("sumr_blank_001", "sumr_blank_048", "sumr_blank_052", 
          "sumr_qc_004", "sumr_qc_015", "sumr_hek_079", "sumr_hek_161",
          "sumr_hek_176", "sumr_hep_006", "sumr_hep_031", "sumr_hep_032"),]
meta$filtered_out <- meta$filtered_out == TRUE

meta <- dplyr::select(meta,-filename)

###
# Extraction Options ( magic, MzT, load )
###
## Load a prior extraction to save time 
#load("mzL.Rdata")
#Magic
mzL <- mzml_extract_magic("data")
mz_obj <- mzL$neg
#save(mzL, file = 'mzL.Rdata')
#load("mzL.Rdata")

## MzT extraction without Magic (Good for debugging).
files <- list.files(path=paste0(getwd(),"/data"), pattern="*.mzML")
files <- paste0("data/",files)
mzT_raw <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
mzT_filt <- pbapply::pblapply(mzT_raw, function(x) mzT_filtering( x, log_name = "manual_neg"))
mzT <- pbapply::pblapply(mzT_filt, function(x) mzT_squash_time(x))
names <- unlist(lapply(files, function(x) strsplit(basename(x),"\\.")[[1]][1]))
for ( x in 1:length(names)) { colnames(mzT[[x]])[2] <- paste0("neg_", names[x])}
#Note, for smaller datasets this is less memory intense than the duckdb join. But be warned, it does not scale well.
mz_raw <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, by= "mz", ...),
    mzT
)
mz_obj <- mzT_binning(mz_raw)

     
# Quality, filtering & post processing
gc()
tech_quality <- mz_quality_magic(mz_obj, meta)
mz_obj <- mz_subtraction(mztools_filter(mz_obj, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_obj, meta,"IS_Blank","type",F,T))
mz_obj <- mz_mass_defect(mz_obj, plot = TRUE)
mz_obj <- mz_filter_magic(mz_obj, min_intensity = 5000, missingness_threshold = 2)
pp_out <- mz_post_magic(mz_obj, metadata = meta, plot = TRUE)

mzlog_analysis_pca(mztools_filter(pp_out$mzLog,meta,"cell","type"), meta)


mzlog_analysis_volcano_magic(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                             mztools_filter(pp_out$mzLog,meta,"HEK293T"))

welch <- mzlog_analysis_welch( mztools_filter(pp_out$mzLog,meta,"HepG2"),
                                 mztools_filter(pp_out$mzLog,meta,"HEK293T"))

fold <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                              mztools_filter(pp_out$mzLog,meta,"HEK293T"))

plot_volcano(welch, fold)

#simpler heatmap for clarity
#mzlog_analysis_heatmap(pp_out$mzLog[1:100,], annotation = NULL, plot_label_size = c(5,2), plot_labels = c(T,T))
mzlog_analysis_heatmap(pp_out$mzLog, metadata = filter(meta, type != "IS_Blank"), plot_label_size = c(5,2), plot_labels = c(T,T))

##RandomForest:Magic
rf_data <- mztools_filter(pp_out$mzLog, meta, c("HepG2","HEK293T"))
#rf_mag <- mzlog_rf_magic(rf_data, meta, "pos")

##RandomForest: without magic
gc()
#High correlation cuttoff due to small sample set
rf_cor <- mzlog_rf_correlation(rf_data, correlation_cutoff = .99)
rf_sel <- mzlog_rf_select(rf_cor, meta)
#RF cannot do much with highly correlated data; needs more samples
rf_obj <- mzlog_rf(rf_data, meta, rfe_obj = rf_sel)
rf_validate <- rf_validate(rf_obj)
rf_permuted <- rf_permuted(rf_obj, "neg")

gc()
anova <- mzlong_analysis_anova(pp_out$mzLong, meta, phenotypes=c("HepG2","HEK293T"), cores=1)

hmdb_rf <- identify_hmdb(rf_mag$mz, c("M-H","M+Na"))
hmdb_anova <- identify_hmdb(anova$imp_mz)
lipids_rf <- identify_lipids(rf_mag$mz)
lipids_rf <- identify_lipids(anova$imp_mz)



