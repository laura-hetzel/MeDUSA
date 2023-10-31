# *** RandomForest Correlation -----------------------------------------------------
#' Find Correlation data within a mzLog_obj
#'
#' @returns Transposed mzlog_obj
#'
#' @export
mzlog_rf_correlation <- function(input_mzlog_obj, correlation_cutoff = 0.75){
  data <- data.frame(t(dplyr::select(input_mzlog_obj,-mz)))
  high_cor_data <- caret::findCorrelation(cor(data), cutoff = correlation_cutoff)
  print(paste("INFO:mzlog_prep: ", length(high_cor_data)," of ", ncol(data) ," MZs are highly correlated ", sep=""))
  data <- data.frame(data[,-high_cor_data])
  colnames(data) <- data["mz",]
  data$sample_name <- rownames(data)
  return(data[-1,])
}

# *** RandomForest Select -----------------------------------------------------
#' Find Correlation data within a mzLog_obj
#'
#' @param correlation_data \cr
#'   DataFrame : from rf_correlation
#' @param metadata \cr
#'   DataFrame: metadata object
#' @param attribute \cr
#'   String: which metadata attribute to compare
#' @param feat_size_seq \cr
#'   Sequence to find optimal "number_of_variables"
#'
#' @returns caret::rfe object
#'
#' @export
rf_select <- function(correlation_data, metadata, attribute = "phenotype", feat_size_seq = seq(50,1000, by=50)) {
  data <- dplyr::left_join(correlation_data, metadata[c("sample_name", attribute)])
  data[[attribute]] <- as.factor(data[[attribute]])

  control <- caret::rfeControl( functions = caret::rfFuncs,
                                method = "repeatedcv",
                                repeats = 5,
                                number = 10)

  feat_select <- caret::rfe(data %>%
                            dplyr::select(-attribute, -"sample_name"),
                            data[[attribute]],
                            rfeControl = control,
                            sizes = seq_size )

  ggplot(data = feat_select, metric = "Accuracy") + theme_bw()
  local.save_plot(paste("RandomForest Accuracy",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
  ggplot(data = feat_select, metric = "Kappa") + theme_bw()
  local.save_plot(paste("RandomForest Kappa",local.mz_polarity_guesser(input_mzlog_obj),sep="-"))
  feat_select
}

# *** RandomForest Train -----------------------------------------------------
#' Train you model data within a mzLog_obj
#'
#' @param mzlog_obj \cr
#'   DataFrame : mzlog
#' @param rfe_obj \cr
#'   carat::rfe_obj (see rf_select)
#' @param metadata \cr
#'   DataFrame: metadata object
#' @param feat_size_seq \cr
#'   Sequence to find optimal "number_of_variables"
#' @param attribute \cr
#'   String: which metadata attribute to compare
#' @param seeds \cr
#'   List: which seeds to use. Also how many runs to do
#' @param cores \cr
#'   Int: can I haz multithreading
#'
#'
#' @export
mzlog_rf_train <- function(mzlog_obj,  metadata, attribute = "phenotype", trees = NULL, rfe_obj = NULL,
                          seeds = c(42,666,314159,1.05457, 998001), ratio = 0.8, cores = 4){
  if (! is.null(rfe_obj)){
    mzLog_obj <- mzLog_obj[ mzLog_obj$mz %in% as.numeric(caret::predictors(rfe_obj)), ]
  }
  if (is.null(trees)){
    if( is.null(rfe_obj)){
      stop("ERROR: mzlog_rf_train: both rfe_obj & trees are null")
    }
    trees <- rfe_obj$bestSubset
  }
  rownames(mzlog_obj) <- mzlog_obj$mz
  data_t <- data.frame(t(dplyr::select(mzlog_obj,-mz)))
  data_t <- tibble::rownames_to_column(data_t,"sample_name")
  data_t <- dplyr::left_join(data_t, metadata[c("sample_name", attribute)])
  data_t <- dplyr::select(data_t,-"sample_name")
  data_t[[attribute]] <- as.factor(data_t[[attribute]])
  
  #Train settings
  mtry <- c(sqrt(ncol(data_t)))
  tunegrid <- expand.grid(.mtry=mtry)
  control <- trainControl(method ='repeatedcv', 
                          number = 10, 
                          repeats = 4, 
                          search = 'grid',
                          allowParallel = TRUE)
  
  cl <- local.export_thread_env(cores, environment(mzlog_rf_train))
  tryCatch({
    out <- pbapply::pblapply(seeds, cl=cl, .run_train)
    out <- do.call(rbind,out)
  }, finally={
    local.kill_threads(cl)
  })
  out
}


.run_train <- function(seed, data_t, ratio, trees, attribute, mtry, tunegrid, control){
  
  #Split
  set.seed(seed)
  split <- caTools::sample.split(data_t$phenotype, SplitRatio = ratio, trees)
  train <- subset(data_t, split == TRUE)
  test  <- subset(data_t, split == FALSE)  
  
  #Train
  
  rf_fit <- train( as.factor(attribute) ~.,
                   data = data_t,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = control,
                   ntree = trees ,
                   na.action = na.exclude)
  
  rf_pred <- predict(rf_fit, test)
  caret::confusionMatrix(rf_pred, as.factor(test_df[[attribute]]))
  
  imp <- caret::varImp(rf_fit)
  imp <- imp$importance
  imp <- rownames_to_column(imp, "mz")
  imp$mz <- as.numeric(gsub("`", "", imp$mz))
  imp <- imp[order(-imp$Overall),]
  imp <- head(imp, n = 100)
}
## set the phenotype to phenotype + a number so that it is unique and can be a row name
#train_neg$phenotype <- paste(train_neg$phenotype, 1:108, sep = "_")
#rownames(train_neg) <- NULL
#train_neg <- column_to_rownames(train_neg, "phenotype")
#
## transpose and add an mz column
#train_neg_t <- as.data.frame(t(train_neg))
#train_neg_t <- rownames_to_column(train_neg_t, "mz")
#imp_var_neg <- train_neg_t %>%
#  filter_all(any_vars(.%in% imp_neg$mz))
#imp_var_neg <- column_to_rownames(imp_var_neg, "mz")
#imp_var_neg_t <- as.data.frame(t(imp_var_neg))
#imp_var_neg_t <- rownames_to_column(imp_var_neg_t, "phenotype")
#imp_var_neg_t$phenotype <- as.factor(ifelse(grepl(imp_var_neg_t$phenotype,
#                                                  pattern = "HepG2"),
#                                            "HepG2", "HEK293T"))
#train_neg <- imp_var_neg_t
## 286 peaks
#
## set the phenotype to phenotype + a number so that it is unique and can be a row name
#train_pos$phenotype <- paste(train_pos$phenotype, 1:66, sep = "_")
#rownames(train_pos) <- NULL
#train_pos <- column_to_rownames(train_pos, "phenotype")
## transpose and add an mz column
#train_pos_t <- as.data.frame(t(train_pos))
#train_pos_t <- rownames_to_column(train_pos_t, "mz")
#
#imp_var_pos <- train_pos_t %>%
#  filter_all(any_vars(.%in% imp_pos$mz))
#imp_var_pos <- column_to_rownames(imp_var_pos, "mz")
#imp_var_pos_t <- as.data.frame(t(imp_var_pos))
#imp_var_pos_t <- rownames_to_column(imp_var_pos_t, "phenotype")
#imp_var_pos_t$phenotype <- as.factor(ifelse(grepl(imp_var_pos_t$phenotype,
#                                                  pattern = "HepG2"),
#                                            "HepG2", "HEK293T"))
#train_pos <- imp_var_pos_t
## 285 peaks
#
#test_neg$phenotype <- paste(test_neg$phenotype, 1:27, sep = "_")
#rownames(test_neg) <- NULL
#test_neg <- column_to_rownames(test_neg, "phenotype")
## transpose and add an mz column
#test_neg_t <- as.data.frame(t(test_neg))
#test_neg_t <- rownames_to_column(test_neg_t, "mz")
#
#test_neg_t <- test_neg_t %>%
#  filter_all(any_vars(.%in% imp_neg$mz))
#test_neg_t <- column_to_rownames(test_neg_t, "mz")
#test_neg <- as.data.frame(t(test_neg_t))
#test_neg <- rownames_to_column(test_neg, "phenotype")
#test_neg$phenotype <- as.factor(ifelse(grepl(test_neg$phenotype,
#                                             pattern = "HepG2"),
#                                       "HepG2", "HEK293T"))
#
#colnames(train_neg) <- as.character(colnames(train_neg))
#
#test_pos$phenotype <- paste(test_pos$phenotype, 1:27, sep = "_")
#rownames(test_pos) <- NULL
#test_pos <- column_to_rownames(test_pos, "phenotype")
## transpose and add an mz column
#test_pos_t <- as.data.frame(t(test_pos))
#test_pos_t <- rownames_to_column(test_pos_t, "mz")
#
#test_pos_t <- test_pos_t %>%
#  filter_all(any_vars(.%in% imp_pos$mz))
#test_pos_t <- column_to_rownames(test_pos_t, "mz")
#test_pos <- as.data.frame(t(test_pos_t))
#test_pos <- rownames_to_column(test_pos, "phenotype")
#test_pos$phenotype <- as.factor(ifelse(grepl(test_pos$phenotype,
#                                             pattern = "HepG2"),
#                                       "HepG2", "HEK293T"))
#
#colnames(train_pos) <- as.character(colnames(train_pos))
#
#set.seed(2017)
#model_neg <- randomForest(x = train_neg[-1],
#                          y = train_neg$phenotype,
#                          data = train_neg,
#                          ntree = 700,
#                          mtry = 20,
#                          importance = TRUE,
#                          proximity = TRUE)
#
#print(model_neg)
#plot(model_neg)
#
#model_pos <- randomForest(x = train_pos[-1],
#                          y = train_pos$phenotype,
#                          data = train_pos,
#                          ntree = 900,
#                          mtry = 20,
#                          importance = TRUE,
#                          proximity = TRUE)
#
#print(model_pos)
#plot(model_pos)
#
## confusion matrix of model - training set
#pred_train_neg <- predict(model_neg, newdata = train_neg)
#confusionMatrix(data = pred_train_neg, reference = train_neg$phenotype)
#pred_train_pos <- predict(model_pos, newdata = train_pos)
#confusionMatrix(data = pred_train_pos, reference = train_pos$phenotype)
#
## prediction of model - test set
#test_neg$phenotype <- as.factor(test_neg$phenotype)
#predict_train_neg <- predict(model_neg, newdata = test_neg, type = "class")
#confusionMatrix(table(data = predict_train_neg, reference = test_neg$phenotype))
#results_neg <- as.data.frame(cbind("actual" = test_neg$phenotype,
#                                   "prediction" = predict_train_neg))
#
#test_pos$phenotype <- as.factor(test_pos$phenotype)
#predict_train_pos <- predict(model_pos, newdata = test_pos, type = "class")
#confusionMatrix(table(data = predict_train_pos, reference = test_pos$phenotype))
#results_pos <- as.data.frame(cbind("actual" = test_pos$phenotype,
#                                   "prediction" = predict_train_pos))
#
#modellist_neg <- data.frame(mtry = 1:200)
#for (i in c(1:200)){
#  set.seed(1230)
#  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
#                      y = train_neg$phenotype,
#                      data = train_neg,
#                      ntree = 500,
#                      mtry = i,
#                      importance = TRUE)
#
#  modellist_neg$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
#}
#
## mtry = 17  has the lowest error rate of approx 5.5%
#
#modellist_pos <- data.frame(mtry = 1:240)
#for (i in c(1:240)){
#  set.seed(1230)
#  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
#                      y = train_neg$phenotype,
#                      data = train_neg,
#                      ntree = 900,
#                      mtry = i,
#                      importance = TRUE)
#
#  modellist_pos$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
#}
## mtry = 3 has the lowest error rate,  approx 5.5%
#
## permutation of data
#permuted_neg <- rbind(train_neg, test_neg)
#x <- c("HepG2", "HEK293T")
## get all permutations
#set.seed(1230)
#false_samples_neg <- sample(x, 135, replace = TRUE)
#
#permuted_neg$phenotype <- as.factor(false_samples_neg)
#permuted_predict_neg <- predict(model_neg, newdata = permuted_neg[-1])
#confusionMatrix(permuted_predict_neg, permuted_neg$phenotype)
## accuracy 0.5111
#results_permuted_neg <- as.data.frame(cbind("actual" = permuted_neg$phenotype,
#                                            "prediction" = permuted_predict_neg))
#roc_permuted_neg <- roc(results_permuted_neg$actual, results_permuted_neg$prediction)
#auc_permuted_neg <- auc(roc_permuted_neg)
#
#permuted_pos <- rbind(train_pos, test_pos)
#x <- c("HepG2", "HEK293T")
#set.seed(1230)
#false_samples_pos <- sample(x, 135, replace = TRUE)
#
#permuted_pos$phenotype <- as.factor(false_samples_pos)
#permuted_predict_pos <- predict(model_pos, newdata = permuted_pos[-1])
#confusionMatrix(permuted_predict_pos, permuted_pos$phenotype)
## accuracy 0.4593
#results_permuted_pos <- as.data.frame(cbind("actual" = permuted_pos$phenotype,
#                                            "prediction" = permuted_predict_pos))
#roc_permuted_pos <- roc(results_permuted_pos$actual, results_permuted_pos$prediction)
#auc_permuted_pos <- auc(roc_permuted_pos)
#
## AUC and ROC, first prep results from above models
#test_neg1$phenotype <- as.factor(ifelse(grepl(test_neg1$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_neg2$phenotype <- as.factor(ifelse(grepl(test_neg2$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_neg3$phenotype <- as.factor(ifelse(grepl(test_neg3$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_neg4$phenotype <- as.factor(ifelse(grepl(test_neg4$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_neg5$phenotype <- as.factor(ifelse(grepl(test_neg5$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#
#test_pos1$phenotype <- as.factor(ifelse(grepl(test_pos1$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_pos2$phenotype <- as.factor(ifelse(grepl(test_pos2$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_pos3$phenotype <- as.factor(ifelse(grepl(test_pos3$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_pos4$phenotype <- as.factor(ifelse(grepl(test_pos4$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#test_pos5$phenotype <- as.factor(ifelse(grepl(test_pos5$phenotype,
#                                              pattern = "HepG2"),
#                                        "HepG2", "HEK293T"))
#
#results_neg1 <- as.data.frame(cbind("actual" = test_neg1$phenotype,
#                                    "prediction" = rf_pred_neg1))
#results_neg2 <- as.data.frame(cbind("actual" = test_neg2$phenotype,
#                                    "prediction" = rf_pred_neg2))
#results_neg3 <- as.data.frame(cbind("actual" = test_neg3$phenotype,
#                                    "prediction" = rf_pred_neg3))
#results_neg4 <- as.data.frame(cbind("actual" = test_neg4$phenotype,
#                                    "prediction" = rf_pred_neg4))
#results_neg5 <- as.data.frame(cbind("actual" = test_neg5$phenotype,
#                                    "prediction" = rf_pred_neg5))
#
#results_pos1 <- as.data.frame(cbind("actual" = test_pos1$phenotype,
#                                    "prediction" = rf_pred_pos1))
#results_pos2 <- as.data.frame(cbind("actual" = test_pos2$phenotype,
#                                    "prediction" = rf_pred_pos2))
#results_pos3 <- as.data.frame(cbind("actual" = test_pos3$phenotype,
#                                    "prediction" = rf_pred_pos3))
#results_pos4 <- as.data.frame(cbind("actual" = test_pos4$phenotype,
#                                    "prediction" = rf_pred_pos4))
#results_pos5 <- as.data.frame(cbind("actual" = test_pos5$phenotype,
#                                    "prediction" = rf_pred_pos5))
#
#roc_neg1 <- roc(results_neg1$actual, results_neg1$prediction)
#auc_neg1 <- auc(roc_neg1)
#
#roc_neg2 <- roc(results_neg2$actual, results_neg2$prediction)
#auc_neg2 <- auc(roc_neg2)
#
#roc_neg3 <- roc(results_neg3$actual, results_neg3$prediction)
#auc_neg3 <- auc(roc_neg3)
#
#roc_neg4 <- roc(results_neg4$actual, results_neg4$prediction)
#auc_neg4 <- auc(roc_neg4)
#
#roc_neg5 <- roc(results_neg5$actual, results_neg5$prediction)
#auc_neg5 <- auc(roc_neg5)
#
#roc_neg_final <- roc(results_neg$actual, results_neg$prediction)
#auc_neg_final <- auc(roc_neg_final)
#
#roc_pos1 <- roc(results_pos1$actual, results_pos1$prediction)
#auc_pos1 <- auc(roc_pos1)
#
#roc_pos2 <- roc(results_pos2$actual, results_pos2$prediction)
#auc_pos2 <- auc(roc_pos2)
#
#roc_pos3 <- roc(results_pos3$actual, results_pos3$prediction)
#auc_pos3 <- auc(roc_pos3)
#
#roc_pos4 <- roc(results_pos4$actual, results_pos4$prediction)
#auc_pos4 <- auc(roc_pos4)
#
#roc_pos5 <- roc(results_pos5$actual, results_pos5$prediction)
#auc_pos5 <- auc(roc_pos5)
#
#roc_pos_final <- roc(results_pos$actual, results_pos$prediction)
#auc_pos_final <- auc(roc_pos_final)
#
## ROC curve of each fold/section
#plot(roc_neg1, col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_neg1, 3))))
#plot(roc_neg2, col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_neg2, 3))))
#plot(roc_neg3, col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_neg3, 3))))
#plot(roc_neg4, col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_neg4, 3))))
#plot(roc_neg5, col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_neg5, 3))))
#plot(roc_neg_final, col = "Blue", main = paste("Final Model(neg), AUC:",
#                                               as.character(round(auc_neg_final, 3))))
#plot(roc_permuted_neg, col = "Green", main = paste("Permuted model(neg), AUC:",
#                                                   as.character(round(auc_permuted_neg, 3))))
#
#plot(roc_pos1, col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_pos1, 3))))
#plot(roc_pos2, col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_pos2, 3))))
#plot(roc_pos3, col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_pos3, 3))))
#plot(roc_pos4, col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_pos4, 3))))
#plot(roc_pos5, col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_pos5, 3))))
#plot(roc_pos_final, col = "Blue", main = paste("Final Model(pos), AUC:",
#                                               as.character(round(auc_pos_final, 3))))
#plot(roc_permuted_pos, col = "Green", main = paste("Permuted model(pos), AUC:",
#                                                   as.character(round(auc_permuted_pos, 3))))
#
## plot the important variables
#varImpPlot(model_neg, col = "Blue", pch = 2, main = "Important neg Variables")
#varImpPlot(model_pos, col = "Blue", pch = 2, main = "Important pos Variables")
#
## MDS plot
#distance_matrix_neg <- as.dist(1 - model_neg$proximity)
#mds_neg <- cmdscale(distance_matrix_neg, eig = TRUE, x.ret = TRUE)
## calculate the percentage of variation that each MDS axis accounts for
#mds_var_neg <- round(mds_neg$eig/sum(mds_neg$eig)*100, 1)
#
#mds_neg_values <- mds_neg$points
#mds_neg_data <- data.frame(Sample = rownames(mds_neg_values),
#                           x = mds_neg_values[,1],
#                           y = mds_neg_values[,2],
#                           Phenotype = train_neg$phenotype)
#
#distance_matrix_pos <- as.dist(1 - model_pos$proximity)
#mds_pos <- cmdscale(distance_matrix_pos, eig = TRUE, x.ret = TRUE)
#mds_var_pos <- round(mds_pos$eig/sum(mds_pos$eig)*100, 1)
#
#mds_pos_values <- mds_pos$points
#mds_pos_data <- data.frame(Sample = rownames(mds_pos_values),
#                           x = mds_pos_values[,1],
#                           y = mds_pos_values[,2],
#                           Phenotype = train_pos$phenotype)
#
#ggplot(mds_neg_data, aes(x = x, y = y, label = Sample)) +
#  geom_point(aes(color = Phenotype), size = 3) +
#  theme_classic() +
#  theme(legend.position = "bottom",) +
#  theme(text = element_text(size = 18)) +
#  xlab(paste("MDS1 - ", mds_var_neg[1], "%", sep = "")) +
#  ylab(paste("MDS2 - ", mds_var_neg[2], "%", sep = "")) +
#  scale_color_manual(breaks = c("HepG2", "HEK293T"),
#                     values = c("red4", "deepskyblue2")) +
#  ggtitle("MDS plot (neg) using random forest proximities")
#
#ggplot(mds_pos_data, aes(x = x, y = y, label = Sample)) +
#  geom_point(aes(color = Phenotype), size = 3) +
#  theme_classic() +
#  theme(legend.position = "bottom",) +
#  theme(text = element_text(size = 18)) +
#  xlab(paste("MDS1 - ", mds_var_pos[1], "%", sep = "")) +
#  ylab(paste("MDS2 - ", mds_var_pos[2], "%", sep = "")) +
#  scale_color_manual(breaks = c("HepG2", "HEK293T"),
#                     values = c("red4", "deepskyblue2")) +
#  ggtitle("MDS plot (pos) using random forest proximities")
#
#save(mds_neg_data, file = "mds_neg_data.R")
#save(mds_pos_data, file = "mds_pos_data.R")
