library('MeDUSA')

# Assumes Docker environment
# set the working directory to where the mzml files are located
setwd("~/local/examples/mzml")
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
## MzT without Magic (Good for debugging)
#files <- list.files(path=getwd("data"), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))

# Quality, filtering & post processing
mz_obj_p <- mzL$pos
tech_quality_p <- mz_quality_magic(mz_obj_p, meta)
mz_obj_p <- mz_subtraction(mztools_filter(mz_obj_p, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_obj_p, meta,"IS_Blank","type",F,T))
mz_obj_p <- mz_mass_defect(mz_obj_p, plot = TRUE)
mz_obj_p <- mz_filter_magic(mz_obj_p, blacklist= F)
pp_out <- mz_post_magic(mz_obj_p, metadata = meta, plot = TRUE)

#mmetadata <- local.meta_polarity_fixer(input_mzlog_obj,metadata)
rownames(input_mzlog_obj) <- input_mzlog_obj$mz
rownames(metadata) <- NULL
t_mz_obj <- scale(t(dplyr::select(input_mzlog_obj,-mz)))
t_mz_obj <- t_mz_obj[!row.names(t_mz_obj) %in% sample_blacklist ,]
t_mz_obj <- merge(x = t_mz_obj,
                  y = subset(tibble::column_to_rownames(metadata, "sample_name"), select =  "phenotype"),
                  by = 0, all.x = TRUE)
rownames(t_mz_obj) <- paste(t_mz_obj$Row.names, t_mz_obj$phenotype, sep = "&")

t_mz_obj <- prcomp(subset(t_mz_obj, select = -c(Row.names, phenotype)))
summary(t_mz_obj)

var_explained <- t_mz_obj$sdev^2/sum(t_mz_obj$sdev^2)

mzlog_analysis_pca(mztools_filter(pp_out$mzLog,meta,"cell","type"), meta)


mzlog_analysis_volcano_magic(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                             mztools_filter(pp_out$mzLog,meta,"HEK293T"))

welch_p <- mzlog_analysis_welch( mztools_filter(pp_out$mzLog,meta,"HepG2"),
                                 mztools_filter(pp_out$mzLog,meta,"HEK293T"))

fold_p <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                              mztools_filter(pp_out$mzLog,meta,"HEK293T"))

plot_volcano(welch_p, fold_p)

##RandomForest:Magic
rf_data <- mztools_filter(pp_out$mzLog, meta, c("HepG2","HEK293T"))
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

