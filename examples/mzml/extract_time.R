library('MeDUSA')

# Assumes Docker environment
# set the working directory to where the mzml files are located
setwd("~/local/examples/mzml")
# introduce the meta data to use for filtering
file_meta <- "meta.xlsx"
meta <- data.frame(readxl::read_excel(file_meta, sheet = "meta", col_names = TRUE))
meta$sample_name <- meta$filename
meta <- dplyr::select(meta,-filename)

files <- list.files(path=paste(getwd(),"data", sep="/"), pattern="*.mzML")

##local.export_thread_env(cores = 4)
cl <- parallel::makeCluster(4, outfile="")
      parallel::clusterExport(cl, varlist = ls(environment()), envir = environment())
      
mzT <- pbapply::pblapply(paste0("data/",files), cl = cl, function(x) MeDUSA::mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
mzT_filtered <- pbapply::pblapply(mzT, cl = cl, function(x) MeDUSA::mzT_filtering( x, log_name = colnames(x)[2]))

parallel::stopCluster(cl)

##extract.mz_format
out <- dplyr::bind_rows(mzT_filtered)
out <- out[order(out$mz),]
out$mz <- round(out$mz, 5)
out <- as.data.frame(out)

out <- MeDUSA::mzT_binning(out)
##Need to export binning method :/

mz_obj_p <- MeDUSA::mz_mass_defect(out, plot = TRUE)
mz_obj_p[is.na(mz_obj_p)] <- 0
mz_obj_p <- MeDUSA::mz_filter_magic(mz_obj_p, missingness_threshold=0.075, min_intensity= 5000, blacklist= F)

mz_obj_pivot <- MeDUSA::mz_post_pivot_longer(mz_obj_p, plot = F)
colnames(mz_obj_pivot) <- c('mz','scanTime','intensity')
mz_obj_pivot$scanTime <- as.numeric(mz_obj_pivot$scanTime)
mz_obj_pivot$intensity <- log2(mz_obj_pivot$intensity)
mz_obj_pivot<- mz_obj_pivot[mz_obj_pivot$intensity != -Inf ,]

subset_mz <- unique(mz_obj_pivot$mz)[seq(1, length(unique(mz_obj_pivot$mz)), 300)]
subset<- mz_obj_pivot[mz_obj_pivot$mz %in% subset_mz,  ]

subsetA<- subset[subset$scanTime < 46.5,  ]
subsetB<- subset[subset$scanTime > 46.5,  ]

ggpubr::ggline(subsetA, "scanTime", "intensity", color = "mz", plot_type = "l",
               numeric.x.axis = T, show.line.label = T)

ggpubr::ggline(subsetB, "scanTime", "intensity", color = "mz", plot_type = "l",
               numeric.x.axis = T, show.line.label = T)
