
#technical_quality_check <- function(dataframe, metadata, mz_name = "mz_s") {
#
#  summary_stats <- data.frame(name = character(length(dataframe)-1 ),
#                              median_mz = numeric(length(dataframe)-1),
#                              min_mz = numeric(length(dataframe)-1),
#                              max_mz = numeric(length(dataframe)-1),
#                              n_peaks = numeric(length(dataframe)-1),
#                              peaks_1k = numeric(length(dataframe)-1),
#                              peaks_10k = numeric(length(dataframe)-1),
#                              peaks_100k = numeric(length(dataframe)-1),
#                              strineqgsAsFactors = F)
#
#
#  # populate the empty data frames
#
#
#  sapply(dataframe[-1], function(column, metadata){
#    column[column>0] %>%
#    dplyr::summarise(measurement = as.numeric(filter(metadata, filename == colnames(column))$measurement),
#              median_mz = median(mz),
#              min_mz = min(mz),
#              max_mz = max(mz),
#              n_peaks = n(),
#              peaks_1k = sum(column > 1000),
#              peaks_10k = sum(column > 10000),
#              peaks_100k = sum(column >100000))
#  })
#
#  for (i in seq_along(dataframe)[-1]) {
#    summary_stats[i-1, ] <- cbind(mz_s = dataframe$mz, dataframe[i]) %>%
#      filter(dataframe[i] > 0) %>%
#      dplyr::summarise(measurement = as.numeric(filter(metadata, filename == colnames(dataframe[i]))$measurement),
#                median_mz = median(mz_s),
#                min_mz = min(mz_s),
#                max_mz = max(mz_s),
#                n_peaks = n(),
#                peaks_1k = sum(dataframe[i] > 1000),
#                peaks_10k = sum(dataframe[i] > 10000),
#                peaks_100k = sum(dataframe[i] >100000))
#  }
#
#  ggplot() +
#    geom_line(data = summary_stats,
#              aes(x = name, y = median_mz, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  #  geom_line(data = summary_stats_p,
#  #            aes(x = name, y = median_mz, colour = "pos", group = 1)) +
#    ggtitle("Median_mz")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = min_mz, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = min_mz, colour = "pos", group = 1)) +
#    ggtitle("Min_mz")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = max_mz, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = max_mz, colour = "pos", group = 1)) +
#    ggtitle("Max_mz")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = n_peaks, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = n_peaks, colour = "pos", group = 1)) +
#    ggtitle("N_peaks")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = peaks_1k, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = peaks_1k, colour = "pos", group = 1)) +
#    ggtitle("Peaks_1k")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = peaks_10k, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = peaks_10k, colour = "pos", group = 1)) +
#    ggtitle("Peaks_10k")
#
#  ggplot() +
#    geom_line(data = summary_stats_n,
#              aes(x = name, y = peaks_100k, colour = "neg", group = 1)) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#    geom_line(data = summary_stats_p,
#              aes(x = name, y = peaks_100k, colour = "pos", group = 1)) +
#    ggtitle("Peaks_100k")
#
#  # remove samples with potentially bad measurements
#  # cell_media_003, sample010, sample006, sample011, and sample012
#  df_pos_clean <- df_pos %>%
#    select(-c("sample011_seg_fullscan_pos",
#              "cell_media_sample003_seg_fullscan_pos",
#              "sample012_seg_fullscan_pos",
#              "cell_media_sample010_seg_fullscan_pos")) %>%
#    filter_at(vars(-matches("mz")), any_vars(. > 0))
#
#  df_pos <- df_pos_clean
#
#  # peaks before tech quality control: 160,857
#  # peaks after tech quality control: 160,837
#  # peaks removed by tech quality control: 20
#  df_neg_clean <- df_neg %>%
#    select(-c("sample011_seg_fullscan_neg",
#              "cell_media_sample003_seg_fullscan_neg",
#              "sample012_seg_fullscan_neg",
#              "cell_media_sample010_seg_fullscan_neg")) %>%
#    filter_at(vars(-matches("mz")), any_vars(. > 0))
#
#  df_neg <- df_neg_clean
#  # peaks before tech quality control: 114,262
#  # peaks after tech quality control: 114,081
#  # peaks removed by tech quality control: 181
#
#  #update df named with removed measurements
#  df_pos_named <- column_to_rownames(df_pos, "mz")
#  df_neg_named <- column_to_rownames(df_neg, "mz")
#
#  # update meta data with removed samples
#  pos_meta$filtered_out[grep("011|003|012|010", pos_meta$filename)] <- "yes"
#  neg_meta$filtered_out[grep("011|003|012|010", neg_meta$filename)] <- "yes"
#
#  # clean up the environment by removing uneccessary data frames
#  rm(df_neg_clean)
#  rm(df_pos_clean)
#}
