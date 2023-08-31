library("sqrlSumr")

##TODO: How to generate this?
#file1_neg <-  "meta_neg.csv"
#neg_meta <- data.frame(read.csv(file1_neg))
#pos_meta <- data.frame(read.csv(file1_pos))

mz_list <- mzml_magic()
mz_pos <- mz_list[['pos']]
mz_neg <- mz_list[['neg']]
