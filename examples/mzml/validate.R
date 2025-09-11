library('MeDUSA')

# Assumes Docker environment
# set the working directory to where the mzml files are located
setwd("~/local/examples/mzml")
# introduce the meta data to use for filtering
file_meta <- "meta.xlsx"
meta <- data.frame(readxl::read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta$sample_name <- meta$filename
meta <- meta[meta$sample_name %in% 
       c("sumr_blank_001", "sumr_blank_048", "sumr_blank_052", 
         "sumr_blank_095", "sumr_qc_004", "sumr_qc_015"),]
meta$filtered_out <- meta$filtered_out == TRUE

meta <- dplyr::select(meta,-filename)

###
# Extraction Options ( magic, MzT, load )
###
## Load a prior extraction to save time 
#load("mzL.Rdata")
#Magic
mzL <- mzml_extract_magic("data")
## MzT extraction without Magic (Good for debugging)
#files <- list.files(path=getwd("data"), pattern="*.mzML")
#mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
#save(mzL, file = 'mzL.Rdata')
#load("mzL.Rdata")
     
# Quality, filtering & post processing
mz_obj_p <- mzL$pos
tech_quality_p <- mz_quality_magic(mz_obj_p, meta)
mz_obj_p <- mz_subtraction(mztools_filter(mz_obj_p, meta,"IS_Blank","type",T,T),
                           mztools_filter(mz_obj_p, meta,"IS_Blank","type",F,T))
mz_obj_p <- mz_mass_defect(mz_obj_p, plot = TRUE)
mz_obj_p <- mz_filter_magic(mz_obj_p, blacklist= F)
pp_out <- mz_post_magic(mz_obj_p, metadata = meta, plot = TRUE)

mzlog_analysis_pca(mztools_filter(pp_out$mzLog,meta,"cell","type"), meta)


mzlog_analysis_volcano_magic(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                             mztools_filter(pp_out$mzLog,meta,"HEK293T"))

welch_p <- mzlog_analysis_welch( mztools_filter(pp_out$mzLog,meta,"HepG2"),
                                 mztools_filter(pp_out$mzLog,meta,"HEK293T"))

fold_p <- mzlog_analysis_fold(mztools_filter(pp_out$mzLog,meta,"HepG2"),
                              mztools_filter(pp_out$mzLog,meta,"HEK293T"))

plot_volcano(welch_p, fold_p)
mzlog_analysis_heatmap(pp_out$mzLog[1:100,], plot_label_size = c(5,2), plot_labels = c(T,T))


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




#### Miscelaneous troubleshooting
library('MeDUSA')
library('duckdb')
setwd("~/local/examples/raw/out")
files <- list.files(getwd(), pattern="*.mzML")
out <- mzml_extract_file(files[1], magic=T, params=list('dbconn'=NULL)) 
mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=T))
mz <- mzml_extract_magic(files)

method <- 'max'
extract.binning(out)

bin <- binning(out$mz, tolerance = 5e-6)
out$mz <- bin
cols <- as.numeric(colnames(select(out, -mz)))
cols <- sprintf('%s("%s")', method, paste(cols ,collapse=paste0('"),',method,'("')))
query <- paste0("SELECT mz, ",cols," FROM out GROUP BY mz" )

setwd("~/local/examples/raw/")
dbconn <- dbConnect(duckdb('tmp.duckdb'))
duckdb::duckdb_register(dbconn, "out", out)
rm(out)

out2 <- dbGetQuery(dbconn, query)

### Probably should delte this; not sure what it's about
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



