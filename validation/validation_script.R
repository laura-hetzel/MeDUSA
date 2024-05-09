
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
library(caret)
library(caTools)
library(randomForest)
library(pROC)

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


# *** 3.4 Random Forest ---------------------------------------------------
# correlate the features
stat_neg_t <- as.data.frame(t(stat_neg))
cor_feat_neg <- cor(stat_neg_t)

# isolate highly correlated, with a correlated ratio of 75%
high_cor_feat_neg <- findCorrelation(cor_feat_neg, cutoff = 0.75)
# 1261 features identified as highly correlated

# identify features to be removed
feat_removal_neg <- as.data.frame(stat_neg_t[, high_cor_feat_neg])
random_neg <- stat_neg_t %>%
  dplyr :: select(-colnames(feat_removal_neg))
# random_neg has 12,221 features
# force the row names to be a column to identify all of the samples
random_neg$samples <- rownames(random_neg)
random_neg <- dplyr:: left_join(random_neg, meta[c("filename", "phenotype")],
                                by = c('samples' = 'filename'))
random_neg$phenotype <- as.factor(random_neg$phenotype)

control_neg <- rfeControl(functions = rfFuncs,
                          method = "repeatedcv",
                          repeats = 5,
                          number = 10)

feat_select_neg <- rfe(random_neg %>%
                         dplyr :: select(-phenotype, -samples),
                       random_neg$phenotype,
                       rfeControl = control_neg,
                       sizes = seq(50,1000, by=50))

ggplot(data = feat_select_neg, metric = "Accuracy") + theme_bw()
ggplot(data = feat_select_neg, metric = "Kappa") + theme_bw()

# mz selection with rows = mz and columns = samples
mz_select_neg <- stat_neg
mz_select_neg <- rownames_to_column(mz_select_neg, "mz")
mz_select_neg <- mz_select_neg %>%
  filter_all(any_vars(.%in% predictors(feat_select_neg)))
# format the data set and add a column to identify the phenotype
mz_select_neg <- column_to_rownames(mz_select_neg, "mz")
mz_select_neg_t <- as.data.frame(t(mz_select_neg))
mz_select_neg_t <- rownames_to_column(mz_select_neg_t, "samples")
rf_neg <- left_join(mz_select_neg_t, select(meta, filename, phenotype),
                    by = c('samples' = 'filename'))
rf_neg <- select(rf_neg, -samples)

# Cross validation

set.seed(42)
split_neg <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg <- subset(rf_neg, split_neg == TRUE)
test_neg <- subset(rf_neg, split_neg == FALSE)

#Train settings
mtry <- c(sqrt(ncol(rf_neg)))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method ='repeatedcv',
                        number = 10,
                        repeats = 4,
                        search = 'grid',
                        allowParallel = TRUE)

rf_fit <- train(as.factor(phenotype) ~.,
                data = train_neg,
                method = 'rf',
                tuneGrid = tunegrid,
                trControl = control,
                ntree = 700,
                na.action = na.exclude)

rf_pred <- predict(rf_fit, test_neg)

confusionMatrix(rf_pred, as.factor(test_neg$phenotype))

# split again for testing, Test number 1
set.seed(111)
split_neg1 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg1 <- subset(rf_neg, split_neg1 == TRUE)
test_neg1 <- subset(rf_neg, split_neg1 == FALSE)

#Train settings same as previous
rf_fit_neg1 <- train(as.factor(phenotype) ~.,
                     data = train_neg1,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 700,
                     na.action = na.exclude)

rf_pred_neg1 <- predict(rf_fit_neg1, test_neg1)
confusionMatrix(rf_pred_neg1, as.factor(test_neg1$phenotype))

# split again for testing, Test number 2
set.seed(2e2)
split_neg2 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg2 <- subset(rf_neg, split_neg2 == TRUE)
test_neg2 <- subset(rf_neg, split_neg2 == FALSE)

#Train settings same as previous

rf_fit_neg2 <- train(as.factor(phenotype) ~.,
                     data = train_neg2,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 700,
                     na.action = na.exclude)

rf_pred_neg2 <- predict(rf_fit_neg2, test_neg2)

confusionMatrix(rf_pred_neg2, as.factor(test_neg2$phenotype))

# split again for testing, Test number 3
set.seed(3e3)
split_neg3 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg3 <- subset(rf_neg, split_neg3 == TRUE)
test_neg3 <- subset(rf_neg, split_neg3 == FALSE)

#Train settings same as previous

#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit_neg3 <- train(as.factor(phenotype) ~.,
                     data = train_neg3,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 700,
                     na.action = na.exclude)

rf_pred_neg3 <- predict(rf_fit_neg3, test_neg3)

confusionMatrix(rf_pred_neg3, as.factor(test_neg3$phenotype))

# split again for testing, Test number 4
set.seed(4e4)
split_neg4 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg4 <- subset(rf_neg, split_neg4 == TRUE)
test_neg4 <- subset(rf_neg, split_neg4 == FALSE)

#Train settings same as previous

#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit_neg4 <- train(as.factor(phenotype) ~.,
                     data = train_neg4,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 700,
                     na.action = na.exclude)

rf_pred_neg4 <- predict(rf_fit_neg4, test_neg4)

confusionMatrix(rf_pred_neg4, as.factor(test_neg4$phenotype))

# split again for testing, Test number 5
set.seed(5e5)
split_neg5 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg5 <- subset(rf_neg, split_neg5 == TRUE)
test_neg5 <- subset(rf_neg, split_neg5 == FALSE)

#Train settings same as previous

#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit_neg5 <- train(as.factor(phenotype) ~.,
                     data = train_neg5,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 700,
                     na.action = na.exclude)

rf_pred_neg5 <- predict(rf_fit_neg5, test_neg5)

confusionMatrix(rf_pred_neg5, as.factor(test_neg5$phenotype))

# Repeat all RFE and random forest for positive mode
# correlate the features
stat_pos_t <- as.data.frame(t(stat_pos))
cor_feat_pos <- cor(stat_pos_t)

# isolate highly correlated, with a correlated ratio of 75%
high_cor_feat_pos <- findCorrelation(cor_feat_pos, cutoff = 0.75)
# 1984 features identified as highly correlated

# identify features to be removed
feat_removal_pos <- as.data.frame(stat_pos_t[, high_cor_feat_pos])
random_pos <- stat_pos_t %>%
  dplyr :: select(-colnames(feat_removal_pos))
# random_pos has 15,637 features
# force the row names to be a column to identify all of the samples
random_pos$samples <- rownames(random_pos)
random_pos <- dplyr:: left_join(random_pos, meta[c("filename", "phenotype")],
                                by = c('samples' = 'filename'))
random_pos$phenotype <- as.factor(random_pos$phenotype)

feat_select_pos <- rfe(random_pos %>%
                         dplyr :: select(-phenotype, -samples),
                       random_pos$phenotype,
                       rfeControl = control_neg,
                       sizes = seq(50,1000, by=50))

ggplot(data = feat_select_pos, metric = "Accuracy") + theme_bw()
ggplot(data = feat_select_pos, metric = "Kappa") + theme_bw()

# mz selection with rows = mz and columns = samples
mz_select_pos <- stat_pos
mz_select_pos <- rownames_to_column(mz_select_pos, "mz")
mz_select_pos <- mz_select_pos %>%
  filter_all(any_vars(.%in% predictors(feat_select_pos)))
# format the data set and add a column to identify the phenotype
mz_select_pos <- column_to_rownames(mz_select_pos, "mz")
mz_select_pos_t <- as.data.frame(t(mz_select_pos))
mz_select_pos_t <- rownames_to_column(mz_select_pos_t, "samples")
rf_pos <- left_join(mz_select_pos_t, select(meta, filename, phenotype),
                    by = c('samples' = 'filename'))
rf_pos <- select(rf_pos, -samples)

# Cross validation

set.seed(42)
split_pos <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos <- subset(rf_pos, split_pos == TRUE)
test_pos <- subset(rf_pos, split_pos == FALSE)

#Train settings
mtry <- c(sqrt(ncol(rf_pos)))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method ='repeatedcv',
                        number = 10,
                        repeats = 4,
                        search = 'grid',
                        allowParallel = TRUE)

rf_fit_pos <- train(as.factor(phenotype) ~.,
                    data = train_pos,
                    method = 'rf',
                    tuneGrid = tunegrid,
                    trControl = control,
                    ntree = 900,
                    na.action = na.exclude)

rf_pred_pos <- predict(rf_fit_pos, test_pos)

confusionMatrix(rf_pred_pos, as.factor(test_pos$phenotype))

# split again for testing, Test number 1
set.seed(111)
split_pos1 <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos1 <- subset(rf_pos, split_pos1 == TRUE)
test_pos1 <- subset(rf_pos, split_pos1 == FALSE)

rf_fit_pos1 <- train(as.factor(phenotype) ~.,
                     data = train_pos1,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 900,
                     na.action = na.exclude)

rf_pred_pos1 <- predict(rf_fit_pos1, test_pos1)

confusionMatrix(rf_pred_pos1, as.factor(test_pos1$phenotype))

# split again for testing, Test number 2
set.seed(2e2)
split_pos2 <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos2 <- subset(rf_pos, split_pos2 == TRUE)
test_pos2 <- subset(rf_pos, split_pos2 == FALSE)

rf_fit_pos2 <- train(as.factor(phenotype) ~.,
                     data = train_pos2,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 900,
                     na.action = na.exclude)

rf_pred_pos2 <- predict(rf_fit_pos2, test_pos2)

confusionMatrix(rf_pred_pos2, as.factor(test_pos2$phenotype))

# split again for testing, Test number 3
set.seed(3e3)
split_pos3 <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos3 <- subset(rf_pos, split_pos3 == TRUE)
test_pos3 <- subset(rf_pos, split_pos3 == FALSE)

rf_fit_pos3 <- train(as.factor(phenotype) ~.,
                     data = train_pos3,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 900,
                     na.action = na.exclude)

rf_pred_pos3 <- predict(rf_fit_pos3, test_pos3)

confusionMatrix(rf_pred_pos3, as.factor(test_pos3$phenotype))

# split again for testing, Test number 4
set.seed(4e4)
split_pos4 <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos4 <- subset(rf_pos, split_pos4 == TRUE)
test_pos4 <- subset(rf_pos, split_pos4 == FALSE)

rf_fit_pos4 <- train(as.factor(phenotype) ~.,
                     data = train_pos4,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 900,
                     na.action = na.exclude)

rf_pred_pos4 <- predict(rf_fit_pos4, test_pos4)

confusionMatrix(rf_pred_pos4, as.factor(test_pos4$phenotype))

# split again for testing, Test number 5
set.seed(5e5)
split_pos5 <- sample.split(rf_pos$phenotype, SplitRatio = 0.8)
train_pos5 <- subset(rf_pos, split_pos5 == TRUE)
test_pos5 <- subset(rf_pos, split_pos5 == FALSE)

rf_fit_pos5 <- train(as.factor(phenotype) ~.,
                     data = train_pos5,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 900,
                     na.action = na.exclude)

rf_pred_pos5 <- predict(rf_fit_pos5, test_pos5)

confusionMatrix(rf_pred_pos5, as.factor(test_pos5$phenotype))

# important variables (features) of each model
# neg
imp_neg1 <- varImp(rf_fit_neg1)
imp_neg1 <- imp_neg1$importance
imp_neg1 <- rownames_to_column(imp_neg1, "mz")
imp_neg1$mz <- as.numeric(gsub("`", "", imp_neg1$mz))
imp_neg1 <- imp_neg1[order(-imp_neg1$Overall),]

imp_neg2 <- varImp(rf_fit_neg2)
imp_neg2 <- imp_neg2$importance
imp_neg2 <- rownames_to_column(imp_neg2, "mz")
imp_neg2$mz <- as.numeric(gsub("`", "", imp_neg2$mz))
imp_neg2 <- imp_neg2[order(-imp_neg2$Overall),]

imp_neg3 <- varImp(rf_fit_neg3)
imp_neg3 <- imp_neg3$importance
imp_neg3 <- rownames_to_column(imp_neg3, "mz")
imp_neg3$mz <- as.numeric(gsub("`", "", imp_neg3$mz))
imp_neg3 <- imp_neg3[order(-imp_neg3$Overall),]

imp_neg4 <- varImp(rf_fit_neg4)
imp_neg4 <- imp_neg4$importance
imp_neg4 <- rownames_to_column(imp_neg4, "mz")
imp_neg4$mz <- as.numeric(gsub("`", "", imp_neg4$mz))
imp_neg4 <- imp_neg4[order(-imp_neg4$Overall),]

imp_neg5 <- varImp(rf_fit_neg5)
imp_neg5 <- imp_neg5$importance
imp_neg5 <- rownames_to_column(imp_neg5, "mz")
imp_neg5$mz <- as.numeric(gsub("`", "", imp_neg5$mz))
imp_neg5 <- imp_neg5[order(-imp_neg5$Overall),]


# isolate the top 100 most important from each model
imp_neg1 <- head(imp_neg1, n = 100)
imp_neg2 <- head(imp_neg2, n = 100)
imp_neg3 <- head(imp_neg3, n = 100)
imp_neg4 <- head(imp_neg4, n = 100)
imp_neg5 <- head(imp_neg5, n = 100)

imp_neg <- rbind(imp_neg1, imp_neg2, imp_neg3, imp_neg4, imp_neg5)

# pos
imp_pos1 <- varImp(rf_fit_pos1)
imp_pos1 <- imp_pos1$importance
imp_pos1 <- rownames_to_column(imp_pos1, "mz")
imp_pos1$mz <- as.numeric(gsub("`", "", imp_pos1$mz))
imp_pos1 <- imp_pos1[order(-imp_pos1$Overall),]

imp_pos2 <- varImp(rf_fit_pos2)
imp_pos2 <- imp_pos2$importance
imp_pos2 <- rownames_to_column(imp_pos2, "mz")
imp_pos2$mz <- as.numeric(gsub("`", "", imp_pos2$mz))
imp_pos2 <- imp_pos2[order(-imp_pos2$Overall),]

imp_pos3 <- varImp(rf_fit_pos3)
imp_pos3 <- imp_pos3$importance
imp_pos3 <- rownames_to_column(imp_pos3, "mz")
imp_pos3$mz <- as.numeric(gsub("`", "", imp_pos3$mz))
imp_pos3 <- imp_pos3[order(-imp_pos3$Overall),]

imp_pos4 <- varImp(rf_fit_pos4)
imp_pos4 <- imp_pos4$importance
imp_pos4 <- rownames_to_column(imp_pos4, "mz")
imp_pos4$mz <- as.numeric(gsub("`", "", imp_pos4$mz))
imp_pos4 <- imp_pos4[order(-imp_pos4$Overall),]

imp_pos5 <- varImp(rf_fit_pos5)
imp_pos5 <- imp_pos5$importance
imp_pos5 <- rownames_to_column(imp_pos5, "mz")
imp_pos5$mz <- as.numeric(gsub("`", "", imp_pos5$mz))
imp_pos5 <- imp_pos5[order(-imp_pos5$Overall),]

# isolate the top 100 most important from each model
imp_pos1 <- head(imp_pos1, n = 100)
imp_pos2 <- head(imp_pos2, n = 100)
imp_pos3 <- head(imp_pos3, n = 100)
imp_pos4 <- head(imp_pos4, n = 100)
imp_pos5 <- head(imp_pos5, n = 100)

imp_pos <- rbind(imp_pos1, imp_pos2, imp_pos3, imp_pos4, imp_pos5)

# set the phenotype to phenotype + a number so that it is unique and can be a row name
train_neg$phenotype <- paste(train_neg$phenotype, 1:108, sep = "_")
rownames(train_neg) <- NULL
train_neg <- column_to_rownames(train_neg, "phenotype")

# transpose and add an mz column
train_neg_t <- as.data.frame(t(train_neg))
train_neg_t <- rownames_to_column(train_neg_t, "mz")
imp_var_neg <- train_neg_t %>%
  filter_all(any_vars(.%in% imp_neg$mz))
imp_var_neg <- column_to_rownames(imp_var_neg, "mz")
imp_var_neg_t <- as.data.frame(t(imp_var_neg))
imp_var_neg_t <- rownames_to_column(imp_var_neg_t, "phenotype")
imp_var_neg_t$phenotype <- as.factor(ifelse(grepl(imp_var_neg_t$phenotype,
                                                  pattern = "HepG2"),
                                            "HepG2", "HEK293T"))
train_neg <- imp_var_neg_t
# 286 peaks

# set the phenotype to phenotype + a number so that it is unique and can be a row name
train_pos$phenotype <- paste(train_pos$phenotype, 1:66, sep = "_")
rownames(train_pos) <- NULL
train_pos <- column_to_rownames(train_pos, "phenotype")
# transpose and add an mz column
train_pos_t <- as.data.frame(t(train_pos))
train_pos_t <- rownames_to_column(train_pos_t, "mz")

imp_var_pos <- train_pos_t %>%
  filter_all(any_vars(.%in% imp_pos$mz))
imp_var_pos <- column_to_rownames(imp_var_pos, "mz")
imp_var_pos_t <- as.data.frame(t(imp_var_pos))
imp_var_pos_t <- rownames_to_column(imp_var_pos_t, "phenotype")
imp_var_pos_t$phenotype <- as.factor(ifelse(grepl(imp_var_pos_t$phenotype,
                                                  pattern = "HepG2"),
                                            "HepG2", "HEK293T"))
train_pos <- imp_var_pos_t
# 285 peaks

test_neg$phenotype <- paste(test_neg$phenotype, 1:27, sep = "_")
rownames(test_neg) <- NULL
test_neg <- column_to_rownames(test_neg, "phenotype")
# transpose and add an mz column
test_neg_t <- as.data.frame(t(test_neg))
test_neg_t <- rownames_to_column(test_neg_t, "mz")

test_neg_t <- test_neg_t %>%
  filter_all(any_vars(.%in% imp_neg$mz))
test_neg_t <- column_to_rownames(test_neg_t, "mz")
test_neg <- as.data.frame(t(test_neg_t))
test_neg <- rownames_to_column(test_neg, "phenotype")
test_neg$phenotype <- as.factor(ifelse(grepl(test_neg$phenotype,
                                             pattern = "HepG2"),
                                       "HepG2", "HEK293T"))

colnames(train_neg) <- as.character(colnames(train_neg))

test_pos$phenotype <- paste(test_pos$phenotype, 1:27, sep = "_")
rownames(test_pos) <- NULL
test_pos <- column_to_rownames(test_pos, "phenotype")
# transpose and add an mz column
test_pos_t <- as.data.frame(t(test_pos))
test_pos_t <- rownames_to_column(test_pos_t, "mz")

test_pos_t <- test_pos_t %>%
  filter_all(any_vars(.%in% imp_pos$mz))
test_pos_t <- column_to_rownames(test_pos_t, "mz")
test_pos <- as.data.frame(t(test_pos_t))
test_pos <- rownames_to_column(test_pos, "phenotype")
test_pos$phenotype <- as.factor(ifelse(grepl(test_pos$phenotype,
                                             pattern = "HepG2"),
                                       "HepG2", "HEK293T"))

colnames(train_pos) <- as.character(colnames(train_pos))

set.seed(2017)
model_neg <- randomForest(x = train_neg[-1],
                          y = train_neg$phenotype,
                          data = train_neg,
                          ntree = 700,
                          mtry = 20,
                          importance = TRUE,
                          proximity = TRUE)

print(model_neg)
plot(model_neg)

model_pos <- randomForest(x = train_pos[-1],
                          y = train_pos$phenotype,
                          data = train_pos,
                          ntree = 900,
                          mtry = 20,
                          importance = TRUE,
                          proximity = TRUE)

print(model_pos)
plot(model_pos)

# confusion matrix of model - training set
pred_train_neg <- predict(model_neg, newdata = train_neg)
confusionMatrix(data = pred_train_neg, reference = train_neg$phenotype)
pred_train_pos <- predict(model_pos, newdata = train_pos)
confusionMatrix(data = pred_train_pos, reference = train_pos$phenotype)

# prediction of model - test set
test_neg$phenotype <- as.factor(test_neg$phenotype)
predict_train_neg <- predict(model_neg, newdata = test_neg, type = "class")
confusionMatrix(table(data = predict_train_neg, reference = test_neg$phenotype))
results_neg <- as.data.frame(cbind("actual" = test_neg$phenotype,
                                   "prediction" = predict_train_neg))

test_pos$phenotype <- as.factor(test_pos$phenotype)
predict_train_pos <- predict(model_pos, newdata = test_pos, type = "class")
confusionMatrix(table(data = predict_train_pos, reference = test_pos$phenotype))
results_pos <- as.data.frame(cbind("actual" = test_pos$phenotype,
                                   "prediction" = predict_train_pos))

modellist_neg <- data.frame(mtry = 1:200)
set.seed(1230)
for (i in c(1:200)){
  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
                      y = train_neg$phenotype,
                      data = train_neg,
                      ntree = 500,
                      mtry = i,
                      importance = TRUE)

  modellist_neg$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}

# mtry = 17  has the lowest error rate of approx 5.5%

modellist_pos <- data.frame(mtry = 1:240)
for (i in c(1:240)){
  set.seed(1230)
  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
                      y = train_neg$phenotype,
                      data = train_neg,
                      ntree = 900,
                      mtry = i,
                      importance = TRUE)

  modellist_pos$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}
# mtry = 3 has the lowest error rate,  approx 5.5%

# permutation of data
permuted_neg <- rbind(train_neg, test_neg)
x <- c("HepG2", "HEK293T")
# get all permutations
set.seed(1230)
false_samples_neg <- sample(x, 135, replace = TRUE)

permuted_neg$phenotype <- as.factor(false_samples_neg)
permuted_predict_neg <- predict(model_neg, newdata = permuted_neg[-1])
confusionMatrix(permuted_predict_neg, permuted_neg$phenotype)
# accuracy 0.5111
results_permuted_neg <- as.data.frame(cbind("actual" = permuted_neg$phenotype,
                                            "prediction" = permuted_predict_neg))
roc_permuted_neg <- roc(results_permuted_neg$actual, results_permuted_neg$prediction)
auc_permuted_neg <- auc(roc_permuted_neg)

permuted_pos <- rbind(train_pos, test_pos)
x <- c("HepG2", "HEK293T")
set.seed(1230)
false_samples_pos <- sample(x, 135, replace = TRUE)

permuted_pos$phenotype <- as.factor(false_samples_pos)
permuted_predict_pos <- predict(model_pos, newdata = permuted_pos[-1])
confusionMatrix(permuted_predict_pos, permuted_pos$phenotype)
# accuracy 0.4593
results_permuted_pos <- as.data.frame(cbind("actual" = permuted_pos$phenotype,
                                            "prediction" = permuted_predict_pos))
roc_permuted_pos <- roc(results_permuted_pos$actual, results_permuted_pos$prediction)
auc_permuted_pos <- auc(roc_permuted_pos)

# AUC and ROC, first prep results from above models
test_neg1$phenotype <- as.factor(ifelse(grepl(test_neg1$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_neg2$phenotype <- as.factor(ifelse(grepl(test_neg2$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_neg3$phenotype <- as.factor(ifelse(grepl(test_neg3$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_neg4$phenotype <- as.factor(ifelse(grepl(test_neg4$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_neg5$phenotype <- as.factor(ifelse(grepl(test_neg5$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))

test_pos1$phenotype <- as.factor(ifelse(grepl(test_pos1$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_pos2$phenotype <- as.factor(ifelse(grepl(test_pos2$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_pos3$phenotype <- as.factor(ifelse(grepl(test_pos3$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_pos4$phenotype <- as.factor(ifelse(grepl(test_pos4$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))
test_pos5$phenotype <- as.factor(ifelse(grepl(test_pos5$phenotype,
                                              pattern = "HepG2"),
                                        "HepG2", "HEK293T"))

results_neg1 <- as.data.frame(cbind("actual" = test_neg1$phenotype,
                                    "prediction" = rf_pred_neg1))
results_neg2 <- as.data.frame(cbind("actual" = test_neg2$phenotype,
                                    "prediction" = rf_pred_neg2))
results_neg3 <- as.data.frame(cbind("actual" = test_neg3$phenotype,
                                    "prediction" = rf_pred_neg3))
results_neg4 <- as.data.frame(cbind("actual" = test_neg4$phenotype,
                                    "prediction" = rf_pred_neg4))
results_neg5 <- as.data.frame(cbind("actual" = test_neg5$phenotype,
                                    "prediction" = rf_pred_neg5))

results_pos1 <- as.data.frame(cbind("actual" = test_pos1$phenotype,
                                    "prediction" = rf_pred_pos1))
results_pos2 <- as.data.frame(cbind("actual" = test_pos2$phenotype,
                                    "prediction" = rf_pred_pos2))
results_pos3 <- as.data.frame(cbind("actual" = test_pos3$phenotype,
                                    "prediction" = rf_pred_pos3))
results_pos4 <- as.data.frame(cbind("actual" = test_pos4$phenotype,
                                    "prediction" = rf_pred_pos4))
results_pos5 <- as.data.frame(cbind("actual" = test_pos5$phenotype,
                                    "prediction" = rf_pred_pos5))

roc_neg1 <- roc(results_neg1$actual, results_neg1$prediction)
auc_neg1 <- auc(roc_neg1)

roc_neg2 <- roc(results_neg2$actual, results_neg2$prediction)
auc_neg2 <- auc(roc_neg2)

roc_neg3 <- roc(results_neg3$actual, results_neg3$prediction)
auc_neg3 <- auc(roc_neg3)

roc_neg4 <- roc(results_neg4$actual, results_neg4$prediction)
auc_neg4 <- auc(roc_neg4)

roc_neg5 <- roc(results_neg5$actual, results_neg5$prediction)
auc_neg5 <- auc(roc_neg5)

roc_neg_final <- roc(results_neg$actual, results_neg$prediction)
auc_neg_final <- auc(roc_neg_final)

roc_pos1 <- roc(results_pos1$actual, results_pos1$prediction)
auc_pos1 <- auc(roc_pos1)

roc_pos2 <- roc(results_pos2$actual, results_pos2$prediction)
auc_pos2 <- auc(roc_pos2)

roc_pos3 <- roc(results_pos3$actual, results_pos3$prediction)
auc_pos3 <- auc(roc_pos3)

roc_pos4 <- roc(results_pos4$actual, results_pos4$prediction)
auc_pos4 <- auc(roc_pos4)

roc_pos5 <- roc(results_pos5$actual, results_pos5$prediction)
auc_pos5 <- auc(roc_pos5)

roc_pos_final <- roc(results_pos$actual, results_pos$prediction)
auc_pos_final <- auc(roc_pos_final)

# ROC curve of each fold/section
plot(roc_neg1, col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_neg1, 3))))
plot(roc_neg2, col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_neg2, 3))))
plot(roc_neg3, col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_neg3, 3))))
plot(roc_neg4, col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_neg4, 3))))
plot(roc_neg5, col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_neg5, 3))))
plot(roc_neg_final, col = "Blue", main = paste("Final Model(neg), AUC:",
                                               as.character(round(auc_neg_final, 3))))
plot(roc_permuted_neg, col = "Green", main = paste("Permuted model(neg), AUC:",
                                                   as.character(round(auc_permuted_neg, 3))))

plot(roc_pos1, col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_pos1, 3))))
plot(roc_pos2, col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_pos2, 3))))
plot(roc_pos3, col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_pos3, 3))))
plot(roc_pos4, col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_pos4, 3))))
plot(roc_pos5, col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_pos5, 3))))
plot(roc_pos_final, col = "Blue", main = paste("Final Model(pos), AUC:",
                                               as.character(round(auc_pos_final, 3))))
plot(roc_permuted_pos, col = "Green", main = paste("Permuted model(pos), AUC:",
                                                   as.character(round(auc_permuted_pos, 3))))

# plot the important variables
varImpPlot(model_neg, col = "Blue", pch = 2, main = "Important neg Variables")
varImpPlot(model_pos, col = "Blue", pch = 2, main = "Important pos Variables")

# MDS plot
distance_matrix_neg <- as.dist(1 - model_neg$proximity)
mds_neg <- cmdscale(distance_matrix_neg, eig = TRUE, x.ret = TRUE)
# calculate the percentage of variation that each MDS axis accounts for
mds_var_neg <- round(mds_neg$eig/sum(mds_neg$eig)*100, 1)

mds_neg_values <- mds_neg$points
mds_neg_data <- data.frame(Sample = rownames(mds_neg_values),
                           x = mds_neg_values[,1],
                           y = mds_neg_values[,2],
                           Phenotype = train_neg$phenotype)

distance_matrix_pos <- as.dist(1 - model_pos$proximity)
mds_pos <- cmdscale(distance_matrix_pos, eig = TRUE, x.ret = TRUE)
mds_var_pos <- round(mds_pos$eig/sum(mds_pos$eig)*100, 1)

mds_pos_values <- mds_pos$points
mds_pos_data <- data.frame(Sample = rownames(mds_pos_values),
                           x = mds_pos_values[,1],
                           y = mds_pos_values[,2],
                           Phenotype = train_pos$phenotype)

ggplot(mds_neg_data, aes(x = x, y = y, label = Sample)) +
  geom_point(aes(color = Phenotype), size = 3) +
  theme_classic() +
  theme(legend.position = "bottom",) +
  theme(text = element_text(size = 18)) +
  xlab(paste("MDS1 - ", mds_var_neg[1], "%", sep = "")) +
  ylab(paste("MDS2 - ", mds_var_neg[2], "%", sep = "")) +
  scale_color_manual(breaks = c("HepG2", "HEK293T"),
                     values = c("red4", "deepskyblue2")) +
  ggtitle("MDS plot (neg) using random forest proximities")

ggplot(mds_pos_data, aes(x = x, y = y, label = Sample)) +
  geom_point(aes(color = Phenotype), size = 3) +
  theme_classic() +
  theme(legend.position = "bottom",) +
  theme(text = element_text(size = 18)) +
  xlab(paste("MDS1 - ", mds_var_pos[1], "%", sep = "")) +
  ylab(paste("MDS2 - ", mds_var_pos[2], "%", sep = "")) +
  scale_color_manual(breaks = c("HepG2", "HEK293T"),
                     values = c("red4", "deepskyblue2")) +
  ggtitle("MDS plot (pos) using random forest proximities")

save(mds_neg_data, file = "mds_neg_data.R")
save(mds_pos_data, file = "mds_pos_data.R")
