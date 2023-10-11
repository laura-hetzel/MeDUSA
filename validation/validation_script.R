
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
  df_pos_named,
  filter(meta, type == "cell" & filtered_out == "no"),
  filter(meta, type == "media" & filtered_out == "no")
)

# peaks before blank subtraction: 116,369
# peaks after blank subtraction: 102,710
# peaks removed with blank subtraction: 13,659

bsub_neg_cells <- blank_subt_local(
  df_neg_named,
  filter(meta, type == "cell" & filtered_out == "no"),
  filter(meta, type == "media" & filtered_out == "no")
)

# peaks before blank subtraction: 121,303
# peaks after blank subtraction: 114,170
# peaks removed with blank subtraction: 7,133


# *** 2.5 Mass Defect filtering -------------------------------------------
filtered_neg_named <- bsub_neg_cells
filtered_pos_named <- bsub_pos_cells

mass_defect_calculation <- function(dataframe) {
  dataframe$MD <-  as.numeric(row.names(dataframe)) %% 1
  dataframe$mz_filter <- 0.00112 * as.numeric(row.names(dataframe)) + 0.01953
  
  # Taken from the else of the SumR thing that was broken anyway 
  #   (revisit if hmdb...McMillan is important)
  md_filtered <- dataframe[which(dataframe$MD <= dataframe$mz_filter), ]
  
  #Ripped from SumR:plot_mz_MD
  mz_removed <- nrow(dataframe) - nrow(md_filtered)
  plot(as.numeric(row.names(md_filtered)), md_filtered$MD,
       cex.axis = 0.8,
       col = alpha("black", 0.5), pch = 20, cex = 0.8,
       ylim = c(0, 1), xlim = c(50, 1200), ylab = "MD", xlab = "m/z",
       main = "Filtered Data", sub = paste("datapoints removed = ", mz_removed),
       cex.lab = 0.8, cex.main = 0.8, cex.sub = 0.8)
  return(select(md_filtered, -MD, -mz_filter))
}

md_filtered_neg_df <- mass_defect_calculation(filtered_neg_named) 
# peaks before md filter: 114,170
# peaks after md filter: 79,284
# peaks removed with md filter: 24,786
md_filtered_pos_df <- mass_defect_calculation(filtered_pos_named) 
# peaks before md filter: 102,710
# peaks after md filter: 74,929
# peaks removed with md filter: 27,781

filtered_neg_named <- md_filtered_neg_df
filtered_pos_named <- md_filtered_pos_df
rm(md_filtered_neg_df)
rm(md_filtered_pos_df)


# *** 2.6 Missingness filter ----------------------------------------------
missingness <- function(dataframe, threshold = 0.1){
  thresh <- threshold * length(dataframe)
  keep_peaks <- dataframe[rowSums( dataframe > 0 ) > thresh,]
}

filtered_neg_named <- missingness(filtered_neg_named)
# peaks before missingness filter: 79,384
# peaks after missingness filter: 30,519
# peaks removed with missingness filter: 48,875

filtered_pos_named <- missingness(filtered_pos_named)
# peaks before missingness filter: 74,929
# peaks after missingness filter: 30,934
# peaks removed with missingness filter: 43,995


# *** 2.7 Low intensity filter --------------------------------------------

neg_int_filter <- filtered_neg_named %>%
  filter_at(vars(-matches("mz")), any_vars( . > 5000))
# peaks before low intensity filter: 30,519
# peaks after low intensity filter: 13,482
# peaks removed by low intensity filter: 17,037

pos_int_filter <- filtered_pos_named %>%
  filter_at(vars(-matches("mz")), any_vars( . > 8000))
# peaks before low intensity filter: 30,934
# peaks after low intensity filter: 17,621
# peaks removed by low intensity filter: 13,313

filtered_neg_named <- neg_int_filter
filtered_pos_named <- pos_int_filter
rm(neg_int_filter)
rm(pos_int_filter)


# *** 2.8 Imputation ------------------------------------------------------
# replace all zeros with a value between 10 and noise
filtered_neg_named[!filtered_neg_named] <- sample(10:5000, 
                                                  sum(!filtered_neg_named),
                                                  replace = TRUE)

filtered_pos_named[!filtered_pos_named] <- sample(10:10000, 
                                                  sum(!filtered_pos_named),
                                                  replace = TRUE)


# *** 2.9 Normalization ---------------------------------------------------
# get the file with the PQN code
source("C:/Users/Laura Hetzel/OneDrive - Universiteit Leiden/Documents/thesis_planning/mimetas_hypoxia/hypoxia_normoxia_mimetas/scripts/quotNorm.R")

# transpose the data frame before normalization
norm_neg <- t(filtered_neg_named)
norm_neg <- quotNorm(norm_neg)

# plot the dilution factor per sample
dilution_neg <- data.frame(norm_neg$dilution)
dilution_neg$sample <- rownames(dilution_neg)
ggplot(dilution_neg, aes(sample, norm_neg.dilution)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90))

dilution_neg_df <- left_join(meta, dilution_neg, by = c('filename' = 'sample'))
ggplot(dilution_neg_df, 
       aes(x = phenotype, y = norm_neg.dilution)) + geom_boxplot()

norm_neg_df <- data.frame(norm_neg$X)
norm_neg_df <- as.data.frame(t(norm_neg_df))

# pos mode
norm_pos <- t(filtered_pos_named)
norm_pos <- quotNorm(norm_pos)

dilution_pos <- data.frame(norm_pos$dilution)
dilution_pos$sample <- rownames(dilution_pos)
ggplot(dilution_pos, aes(sample, norm_pos.dilution)) + geom_point() + theme(axis.text.x = element_text(angle = 90))


dilution_pos_df <- left_join(meta, dilution_pos, by = c('filename' = 'sample'))

ggplot(dilution_pos_df, 
       aes(x = phenotype, y = norm_pos.dilution)) + geom_boxplot()

norm_pos_df <- data.frame(norm_pos$X)
norm_pos_df <- as.data.frame(t(norm_pos_df))

# remove measurements that were flagged by normalization plots and raw data review
# revealed a cell was not measured.
norm_pos_df <- norm_pos_df %>%
  select(-c("sumr_hep_081")) %>%
  filter_at(vars(-matches("mz")), any_vars(. > 0))

norm_neg_df <- norm_neg_df %>%
  select(-c("sumr_hep_081")) %>%
  filter_at(vars(-matches("mz")), any_vars(. > 0))

# update meta data with removed samples
meta$filtered_out[grep("081", meta$filename)] <- "yes"

#rearrange the normalized data sets with pivot longer
norm_neg_df$mz <- row.names(norm_neg_df)
norm_neg_long <- norm_neg_df %>% 
  pivot_longer(!mz, names_to = "sample", values_to = "intensity")

norm_pos_df$mz <- row.names(norm_pos_df)
norm_pos_long <- norm_pos_df %>% 
  pivot_longer(!mz, names_to = "sample", values_to = "intensity")


# *** 2.10 Log Transform --------------------------------------------------
norm_pos_long$log <- log2(norm_pos_long$intensity)
norm_neg_long$log <- log2(norm_neg_long$intensity)

# remove the X from the mz value
norm_pos_long$mz <- gsub("X", "", norm_pos_long$mz)
norm_neg_long$mz <- gsub("X", "", norm_neg_long$mz)

ggplot(norm_neg_long, aes(x = sample, y = log)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Filtered, Normalized, Logged, Negative")

ggplot(norm_pos_long, aes(x = sample, y = log)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Filtered, Normalized, Logged, Positive")


# 3. Statistical Analysis -------------------------------------------------
stat_neg <- subset(norm_neg_df, select = -c(mz))
row.names(stat_neg) <- gsub("X", "", rownames(stat_neg))
stat_neg <- log2(stat_neg)
save(stat_neg, file = "sumr_neg_filtered_for_stats_df.Rdata")
stat_pos <- subset(norm_pos_df, select = -c(mz))
row.names(stat_pos) <- gsub("X", "", rownames(stat_pos))
stat_pos <- log2(stat_pos)
save(stat_pos, file = "stem_pos_filtered_for_stats_df.Rdata")
save(meta, file = "sumr_filtered_meta.Rdata")


# *** 3.1 PCA -------------------------------------------------------------
# neg mode
pca_input_neg <- as.data.frame(t(stat_neg))
pca_input_neg <- scale(pca_input_neg)

pca_input_neg <- merge(x = pca_input_neg, 
                       y = subset(column_to_rownames(meta, "filename"), 
                                  select = "phenotype"),
                       by = 0, all.x = TRUE)
rownames(pca_input_neg) <- paste(pca_input_neg$Row.names, pca_input_neg$phenotype, 
                                 sep = "&")

pca_neg <- prcomp(subset(pca_input_neg, select = -c(Row.names, phenotype)))
summary(pca_neg)

var_explained_n <- pca_neg$sdev^2/sum(pca_neg$sdev^2)

pca_neg$x %>%
  as.data.frame %>%
  rownames_to_column("sample_phenotype") %>%
  separate(sample_phenotype, c("sample", "phenotype"), "&") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = phenotype)) +
  labs(x=paste0("PC1: ",round(var_explained_n[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained_n[2]*100,1),"%")) +
  theme(legend.position="top")

# pos mode
pca_input_pos <- as.data.frame(t(stat_pos))
pca_input_pos <- scale(pca_input_pos)

pca_input_pos <- merge(x = pca_input_pos, 
                       y = subset(column_to_rownames(meta, "filename"), 
                                  select = "phenotype"),
                       by = 0, all.x = TRUE)
rownames(pca_input_pos) <- paste(pca_input_pos$Row.names, pca_input_pos$phenotype, 
                                 sep = "&")

pca_pos <- prcomp(subset(pca_input_pos, select = -c(Row.names, phenotype)))
summary(pca_pos)

var_explained_p <- pca_pos$sdev^2/sum(pca_pos$sdev^2)

pca_pos$x %>%
  as.data.frame %>%
  rownames_to_column("sample_phenotype") %>%
  separate(sample_phenotype, c("sample", "phenotype"), "&") %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = phenotype)) +
  labs(x=paste0("PC1: ",round(var_explained_n[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained_n[2]*100,1),"%")) +
  theme(legend.position="top")


# *** 3.2 T-test ----------------------------------------------------------
welch <- function(df, meta){
  samples_hep <- filter(meta, phenotype == "HepG2" &
                           type == "cell" &
                           filtered_out == "no")
  samples_hek <- filter(meta, phenotype == "HEK293T" &
                            type == "cell" &
                            filtered_out == "no")
  
  a <- data.frame(p    = rep(Inf, nrow(df)), 
                  p_05 = rep(FALSE, nrow(df)),
                  p_10 = rep(FALSE, nrow(df)),
                  p_15 = rep(FALSE, nrow(df)),
                  row.names = row.names(df))
  for (i in 1:nrow(df)){
    # create a new data frame to write the t-test values to
    a$p[i] <- t.test(select(df[i,], samples_hep$filename),
                     select(df[i,], samples_hek$filename))$p.value
  }
  a$p_05 <- a$p < 0.05
  a$p_10 <- a$p < 0.10
  a$p_15 <- a$p < 0.15
  
  a
}

welch_neg <- welch(stat_neg, meta)
welch_pos <- welch(stat_pos, meta)

fdr_neg <- as.data.frame(p.adjust(welch_neg$p), method = p.adjust.methods("fdr"))
# all values are 1, so not helpful
fdr_pos <- as.data.frame(p.adjust(welch_pos$p), method = p.adjust.methods("fdr"))
# all values are 1, so not helpful

p05neg <- select(welch_neg[welch_neg$p_05,],p)
# 492 mz with a p < 0.05
p05pos <- select(welch_pos[welch_pos$p_05,],p)
# 703 mz with a p < 0.05


# *** 3.3 Volcano Plots ---------------------------------------------------

# ***** Fold Change -------------------------------------------------------
fold_change <- function(input_data, num_samples, den_samples) {
  input_data <- as.data.frame(input_data)
  
  num_subset <- input_data[,num_samples$filename]
  den_subset <- input_data[,den_samples$filename]
  input_data$fold <- (apply(num_subset, 1, mean) /
                        apply(den_subset, 1, mean))
  
  df_out <- select(input_data, fold)
}

neg_fold <- fold_change(
  stat_neg,
  filter(meta, phenotype == "HepG2" & type == "cell" & filtered_out == "no"),
  filter(meta, phenotype != "HepG2" & type == "cell"  & filtered_out == "no"))

pos_fold <- fold_change(
  stat_pos,
  filter(meta, phenotype == "HepG2" & type == "cell" & filtered_out == "no"),
  filter(meta, phenotype != "HepG2" & type == "cell"  & filtered_out == "no"))


# ***** fold vs welch volcano ---------------------------------------------

volcano_input_n <- cbind(rownames(welch_neg),welch_neg$p, neg_fold$fold)
colnames(volcano_input_n) <- c("mz", "pvalue", "fold")
volcano_input_n <- data.frame(volcano_input_n)
volcano_input_n$mz <- as.numeric(volcano_input_n$mz)
volcano_input_n$pvalue <- as.numeric(volcano_input_n$pvalue)
volcano_input_n$fold <- as.numeric(volcano_input_n$fold)

volcano_input_n$diff <- "NO"
volcano_input_n$diff[log2(volcano_input_n$fold) > 0.6 & 
                       -log10(volcano_input_n$pvalue) > -log10(0.05)] <- "UP"
volcano_input_n$diff[log2(volcano_input_n$fold) < -0.6 & 
                       -log10(volcano_input_n$pvalue) > -log10(0.05)] <- "DOWN"


ggplot(data = volcano_input_n, aes(x = log2(fold),
                                  y = -log10(pvalue),
                                  col = diff)) +
  geom_point() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  ggtitle("Negative Mode")

volcano_input_p <- cbind(rownames(welch_pos),welch_pos$p, pos_fold$fold)
colnames(volcano_input_p) <- c("mz", "pvalue", "fold")
volcano_input_p <- data.frame(volcano_input_p)
volcano_input_p$mz <- as.numeric(volcano_input_p$mz)
volcano_input_p$pvalue <- as.numeric(volcano_input_p$pvalue)
volcano_input_p$fold <- as.numeric(volcano_input_p$fold)

volcano_input_p$diff <- "NO"
volcano_input_p$diff[log2(volcano_input_p$fold) > 0.6 & 
                       -log10(volcano_input_p$pvalue) > -log10(0.05)] <- "UP"
volcano_input_p$diff[log2(volcano_input_p$fold) < -0.6 & 
                       -log10(volcano_input_p$pvalue) > -log10(0.05)] <- "DOWN"


ggplot(data = volcano_input_p, aes(x = log2(fold),
                                   y = -log10(pvalue),
                                   col = diff)) +
  geom_point() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  ggtitle("Positive Mode")