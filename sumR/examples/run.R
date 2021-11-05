library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(rfUtilities)
library(ElemStatLearn)

### Script
## merge data
options(scipen = 999)
options(digits = 10)
path <- "E:/code/CTC_HU/renamed_rawdata/renamed_pos"
setwd(path)
files <- list.files(path = path)
file_names <- stringr::str_remove(string = files, pattern = ".txt")
file_list <- lapply(files, function(x) read.table(x, col.names = c("mz", "intensity")))

nn <- sapply(file_list, nrow) # each sample's row number
nonEmpty <- nn != 0L # if there's any empty sample; return logic vector
samples <- rep(file_names, nn)
mass <- unname(unlist((lapply(file_list[nonEmpty], function(x) x$mz)), recursive = FALSE, use.names = FALSE))
intensities <- unlist(lapply(file_list[nonEmpty], function(x) x$intensity), recursive = FALSE, use.names = FALSE)
s <- sort.int(mass, index.return = TRUE)
mass <- s$x
intensities <- intensities[s$ix]
samples <- samples[s$ix]

msdata <- data.frame(mz = mass, sample = samples, intensity = intensities)


# withIDdata <- mergeData("E:/code/CTC_HU/renamed_rawdata/HB")
#
# ## transform the data format into the one needed for densityclust()
# datalonger <- withIDdata %>%
#   pivot_longer(!mz, names_to = "sample", values_to = "intensity")
#
# ## remove rows with 0 intensity
# msdata <<- datalonger[(which(datalonger$intensity != 0)), ]

## remove noise
#msdata <- noiseRemove2(msdata, estimatedNoise = 300, SNRlimit = 10)

# initialize a new data frame
sample_num <- length(files)
# AlignedMS <- as.data.frame(t(data.frame(rep(0,sample_num+1))))
# colnames(AlignedMS) <- c("aligned", file_names)

# create bins first
bin <- binning(msdata$mz, tolerance = 2000, times = 4) #two type of conditions can be chosen; tolerance & times belong to different conditions
bin_df <- data.frame(left = bin$left, right = bin$right)
bin_df <- arrange(bin_df, bin_df$left)
bin_split <- split(bin_df, 1:nrow(bin_df))

# check the ppm between the first and last mz of the bin
# bin_df$ppm <- abs(ppm_calc(bin_df$left,bin_df$right))
# bin_df[which(bin_df$ppm<5),]

# processing without parallel
# list_files_aligned <- lapply(bin_split, cluster_align)

# processing with parallel
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl,as.list(unique(ls(environment(cluster_align)))), envir = environment(cluster_align))

list_files_aligned_pl <- parLapply(cl, bin_split, cluster_align)

stopCluster(cl)

# list to dataframe
Alligned_PM <- do.call(rbind, list_files_aligned_pl)
colnames(Alligned_PM) <- c("mz", file_names)





