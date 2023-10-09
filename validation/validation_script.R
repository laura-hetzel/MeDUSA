
# 1. General Info ---------------------------------------------------------
# The purpose of this script is to process and filter data in a way that is 
# standard such that the results may be compared to results obtained with the 
# sumr package. Ultimately, this script will help determine if sumr was written
# in such a way as to conform to industry standards.

# *** 1.1 Sample Info -----------------------------------------------------
# cells - HEPG2 or Hek293T
# blank_solvent - sample of only the ionization solvent to serve as a baseline
#                 for which metabolites are in the solvent itself
# blank_media - media sampled from the either set of cells
# QC - quality control sample, diluted calibration mix

# *** 1.2 Dependencies ----------------------------------------------------
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggbiplot)
library(dplyr)

# 2. Data Processing ------------------------------------------------------


# *** 2.1 Import all files ------------------------------------------------
file0 <- "C:/Users/Laura Hetzel/OneDrive - Universiteit Leiden/Documents/thesis_planning/cellular_signal_r_package/sumr_validation/sumr_validation_master_mz.xlsx"
file1 <- "C:/Users/Laura Hetzel/OneDrive - Universiteit Leiden/Documents/thesis_planning/cellular_signal_r_package/sumr_validation/meta.xlsx"

df_neg <- data.frame(read_excel(file0, sheet = "neg", col_names = TRUE))
# 121,423 peaks

df_pos <- data.frame(read_excel(file0, sheet = "pos", col_names = TRUE))
# 116,407 peaks

colnames(df_neg) <- gsub("\\.txt$", "", colnames(df_neg))
colnames(df_pos) <- gsub("\\.txt$", "", colnames(df_pos))
colnames(df_neg) <- gsub("\\_neg$", "", colnames(df_neg))
colnames(df_pos) <- gsub("\\_pos$", "", colnames(df_pos))

# remove the mz column and instead have the mz be the row names
df_pos_mznamed <- column_to_rownames(df_pos, "mz")
df_neg_mznamed <- column_to_rownames(df_neg, "mz")

# assign the metadata to a data frame
meta <- data.frame(read_excel(file1, sheet = "meta", col_names = TRUE))


# *** 2.2 Technical quality check -----------------------------------------
# establish empty dataframes
summary_stats_p <- data.frame(name = character(length(df_pos)-1 ), 
                              median_mz = numeric(length(df_pos)-1), 
                              min_mz = numeric(length(df_pos)-1), 
                              max_mz = numeric(length(df_pos)-1),
                              n_peaks = numeric(length(df_pos)-1),
                              peaks_1k = numeric(length(df_pos)-1),
                              peaks_10k = numeric(length(df_pos)-1),
                              peaks_100k = numeric(length(df_pos)-1),
                              stringsAsFactors = F)

summary_stats_n <- data.frame(name = character(length(df_neg)-1 ), 
                              median_mz = numeric(length(df_neg)-1), 
                              min_mz = numeric(length(df_neg)-1), 
                              max_mz = numeric(length(df_neg)-1),
                              n_peaks = numeric(length(df_neg)-1),
                              peaks_1k = numeric(length(df_neg)-1),
                              peaks_10k = numeric(length(df_neg)-1),
                              peaks_100k = numeric(length(df_neg)-1),
                              stringsAsFactors = F)

# populate the empty data frames
for (i in seq_along(df_pos)[-1]) {
  summary_stats_p[i-1, ] <- cbind(mz_s = df_pos$mz, df_pos[i]) %>%
    filter(df_pos[i] > 0) %>%
    dplyr::summarise(measurement = as.numeric(filter(meta, filename == colnames(df_pos[i]))$measurement), 
                     median_mz = median(mz_s), 
                     min_mz = min(mz_s),
                     max_mz = max(mz_s), 
                     n_peaks = n(),
                     peaks_1k = sum(df_pos[i] > 1000),
                     peaks_10k = sum(df_pos[i] > 10000),
                     peaks_100k = sum(df_pos[i] >100000))
}
