
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

for (i in seq_along(df_neg)[-1]) {
  summary_stats_n[i-1, ] <- cbind(mz_s = df_neg$mz, df_neg[i]) %>%
    filter(df_neg[i] > 0) %>%
    dplyr::summarise(measurement = as.numeric(filter(meta, filename == colnames(df_neg[i]))$measurement), 
                     median_mz = median(mz_s), 
                     min_mz = min(mz_s),
                     max_mz = max(mz_s), 
                     n_peaks = n(),
                     peaks_1k = sum(df_neg[i] > 1000),
                     peaks_10k = sum(df_neg[i] > 10000),
                     peaks_100k = sum(df_neg[i] >100000))
}

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = median_mz, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = median_mz, colour = "pos", group = 1)) +
  ggtitle("Median_mz")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = min_mz, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = min_mz, colour = "pos", group = 1)) +
  ggtitle("Min_mz")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = max_mz, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = max_mz, colour = "pos", group = 1)) +
  ggtitle("Max_mz")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = n_peaks, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = n_peaks, colour = "pos", group = 1)) +
  ggtitle("N_peaks")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = peaks_1k, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = peaks_1k, colour = "pos", group = 1)) +
  ggtitle("Peaks_1k")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = peaks_10k, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = peaks_10k, colour = "pos", group = 1)) +
  ggtitle("Peaks_10k")

ggplot() +
  geom_line(data = summary_stats_n, 
            aes(x = name, y = peaks_100k, colour = "neg", group = 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_line(data = summary_stats_p, 
            aes(x = name, y = peaks_100k, colour = "pos", group = 1)) +
  ggtitle("Peaks_100k")

# remove samples with potentially bad measurements
# sumr_media_100, sumr_hek_61 both had indications of severe ion suppression
# due to contamination; sumr_media_147 consistently had low signal in pos mode
df_pos_clean <- df_pos %>%
  select(-c("sumr_media_100", 
            "sumr_hek_061",
            "sumr_media_147")) %>%
  filter_at(vars(-matches("mz")), any_vars(. > 0))

df_pos <- df_pos_clean

# peaks before tech quality control: 116,407
# peaks after tech quality control: 116,369
# peaks removed by tech quality control: 38

df_neg_clean <- df_neg %>%
  select(-c("sumr_media_100", 
            "sumr_hek_061",
            "sumr_media_147")) %>%
  filter_at(vars(-matches("mz")), any_vars(. > 0))

df_neg <- df_neg_clean
# peaks before tech quality control: 121,423
# peaks after tech quality control: 121,303
# peaks removed by tech quality control: 120

#update df named with removed measurements
df_pos_named <- column_to_rownames(df_pos, "mz")
df_neg_named <- column_to_rownames(df_neg, "mz")

# update meta data with removed samples
meta$filtered_out[grep("100|061|147", meta$filename)] <- "yes"

# clean up the environment by removing uneccessary data frames
rm(df_neg_clean)
rm(df_pos_clean)


# *** 2.3 Alignment/Binning -----------------------------------------------
# Binning performed in MarkerView software, 5 ppm centroiding option


# *** 2.4 Blank Subtraction -----------------------------------------------
blank_subt_local <- function(input_data, sample_filenames, blank_filenames, blank_thresh = 3) {
  input_data <- as.data.frame(input_data)
  
  blanks_subset <- input_data[,blank_filenames$filename]
  blanks_subset$threshold <- apply(blanks_subset, 1, median) * blank_thresh
  
  samples_subset <- input_data[,sample_filenames$filename]
  
  blank_applyer <- function(input_data, blank = blanks_subset$threshold){
    (input_data > blank) * input_data
  }
  samples_subset <- as.data.frame(sapply(samples_subset, blank_applyer))
  row.names(samples_subset) <- row.names(input_data)
  samples_subset[rowSums(samples_subset) > 0,]
}

bsub_pos_cells <- blank_subt_local(
  df_pos,
  filter(meta, type == "cell" & filtered_out == "no"),
  filter(meta, type == "media" & filtered_out == "no")
)

# peaks before blank subtraction: 116,369
# peaks after blank subtraction: 102,710
# peaks removed with blank subtraction: 13,659

bsub_neg_cells <- blank_subt_local(
  df_neg,
  filter(meta, type == "cell" & filtered_out == "no"),
  filter(meta, type == "media" & filtered_out == "no")
)

# peaks before blank subtraction: 121,303
# peaks after blank subtraction: 114,170
# peaks removed with blank subtraction: 7,133

