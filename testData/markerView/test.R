
library("sqrlSumr")
#source("../../sqrlSumr/R/mz_4_filter.R")


file0_neg <-  "neg.csv"
file1_neg <-  "meta_neg.csv"
mz_obj <- data.frame(read.csv(file0_neg))
neg_meta <- data.frame(read.csv(file1_neg))
colnames(mz_obj) <- gsub("\\.txt$", "", colnames(mz_obj))
mz_obj <- tibble::column_to_rownames(mz_obj, "mz")
neg_meta <- dplyr::rename(neg_meta,c(sample_name=filename))

mz_obj2 <- sqrlSumr::mz_subtraction(mz_obj,
            dplyr::filter(neg_meta, type != "solvent" & filtered_out == "no"),
            dplyr::filter(neg_meta, type == "solvent" & filtered_out == "no"))

mz_obj3 <- sqrlSumr::mz_subtraction(mz_obj2,
            dplyr::filter(neg_meta, type == "cell" & filtered_out == "no"),
            dplyr::filter(neg_meta, type == "media_cell" & filtered_out == "no"))

mz_obj4 <- sqrlSumr::mz_mass_defect(mz_obj3)

#misssingness, LowIntensity
mz_obj5 <- sqrlSumr::mz_filter_magic(mz_obj4)

mz_mag <- sqrlSumr::mz_pp_imputation(mz_obj5)
mz_mag <- sqrlSumr::mz_pp_normalization(mz_mag, neg_meta)
mz_mag <- log2(mz_mag)

#imputation, normaliztion, pivot_longer, log_transform
c(mzl_obj, mz_obj6) <- sqrlSumr::mz_post_process_magic(mz_obj5, neg_meta)

bad_neg_samples <- c("sample155_seg_fullscan_neg",
                     "sample075_seg_fullscan_neg",
                     "sample168_seg_fullscan_neg",
                     "sample181_seg_fullscan_neg",
                     "sample164_seg_fullscan_neg")
sqrlSumr::mz_pca(mz_mag, neg_meta, bad_neg_samples )

samples_diff <- filter(neg_meta, phenotype == "differentiated" &
                         type == "cell" &
                         filtered_out == "no")
samples_naive <- filter(neg_meta, phenotype == "naive" &
                          type == "cell" &
                          filtered_out == "no")

sqrlSumr::mz_welch(mz_mag, neg_meta , samples_diff, samples_naive)
