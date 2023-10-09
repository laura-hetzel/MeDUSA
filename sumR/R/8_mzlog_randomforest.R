# *** 3.5 random forest ---------------------------------------------------

# correlate the features
stat_neg_t <- as.data.frame(t(stat_neg))
cor_feat_neg <- cor(stat_neg_t)
save(cor_feat_neg, file = "stem_neg_correlated_features.Rdata")

# isolate highly correlated, with a correlated ratio of 75%
high_cor_feat_neg <- findCorrelation(cor_feat_neg, cutoff = 0.75)
# 651 features identified as highly correlated

# identify features to be removed
feat_removal_neg <- as.data.frame(stat_neg_t[, high_cor_feat_neg])
random_neg <- stat_neg_t %>%
  dplyr :: select(-colnames(feat_removal_neg))
# random_neg has 6,675 features
# force the row names to be a column to identify all of the samples
random_neg$samples <- rownames(random_neg)
random_neg <- dplyr:: left_join(random_neg, neg_meta[c("filename", "phenotype")],
                                by = c('samples' = 'filename'))
random_neg$phenotype <- as.factor(random_neg$phenotype)

#subsets_cor <- c(150:300, 500)

control_neg <- rfeControl(functions = rfFuncs,
                          method = "repeatedcv",
                          repeats = 5,
                          number = 10)

feat_select_neg <- rfe(random_neg %>%
                         dplyr :: select(-phenotype, -samples),
                       random_neg$phenotype,
                       rfeControl = control_neg,
                       sizes = seq(50,500, by=50))

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
rf_neg <- left_join(mz_select_neg_t, select(neg_meta, filename, phenotype),
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
#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit <- train(as.factor(phenotype) ~.,
                data = train_neg,
                method = 'rf',
                tuneGrid = tunegrid,
                trControl = control,
                ntree = 500,
                na.action = na.exclude)

rf_pred <- predict(rf_fit, test_neg)

confusionMatrix(rf_pred, as.factor(test_neg$phenotype))

# split again for testing, Test number 1
set.seed(111)
split_neg1 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg1 <- subset(rf_neg, split_neg1 == TRUE)
test_neg1 <- subset(rf_neg, split_neg1 == FALSE)

#Train settings same as previous

#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit_neg1 <- train(as.factor(phenotype) ~.,
                data = train_neg1,
                method = 'rf',
                tuneGrid = tunegrid,
                trControl = control,
                ntree = 500,
                na.action = na.exclude)

rf_pred_neg1 <- predict(rf_fit_neg1, test_neg1)

# split again for testing, Test number 2
set.seed(2e2)
split_neg2 <- sample.split(rf_neg$phenotype, SplitRatio = 0.8)
train_neg2 <- subset(rf_neg, split_neg2 == TRUE)
test_neg2 <- subset(rf_neg, split_neg2 == FALSE)

#Train settings same as previous

#rf_fit <- train(as.factor(phenotype) ~., data = test_neg, method= 'rf')
rf_fit_neg2 <- train(as.factor(phenotype) ~.,
                     data = train_neg2,
                     method = 'rf',
                     tuneGrid = tunegrid,
                     trControl = control,
                     ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
                     na.action = na.exclude)

rf_pred_neg5 <- predict(rf_fit_neg5, test_neg5)

confusionMatrix(rf_pred_neg5, as.factor(test_neg5$phenotype))

# Repeat all RFE and random forest for positive mode
# correlate the features
stat_pos_t <- as.data.frame(t(stat_pos))
cor_feat_pos <- cor(stat_pos_t)
save(cor_feat_pos, file = "stem_pos_correlated_features.Rdata")

# isolate highly correlated, with a correlated ratio of 75%
high_cor_feat_pos <- findCorrelation(cor_feat_pos, cutoff = 0.75)
# 644 features identified as highly correlated

# identify features to be removed
feat_removal_pos <- as.data.frame(stat_pos_t[, high_cor_feat_pos])
random_pos <- stat_pos_t %>%
  dplyr :: select(-colnames(feat_removal_pos))
# random_pos has 9,132 features
# force the row names to be a column to identify all of the samples
random_pos$samples <- rownames(random_pos)
random_pos <- dplyr:: left_join(random_pos, pos_meta[c("filename", "phenotype")],
                                by = c('samples' = 'filename'))
random_pos$phenotype <- as.factor(random_pos$phenotype)

feat_select_pos <- rfe(random_pos %>%
                         dplyr :: select(-phenotype, -samples),
                       random_pos$phenotype,
                       rfeControl = control_neg,
                       sizes = seq(50,500, by=50))

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
rf_pos <- left_join(mz_select_pos_t, select(pos_meta, filename, phenotype),
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
                ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
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
                     ntree = 500,
                     na.action = na.exclude)

rf_pred_pos5 <- predict(rf_fit_pos5, test_pos5)

confusionMatrix(rf_pred_pos5, as.factor(test_pos5$phenotype))

# important variables (features) of each model
# neg
imp_neg1 <- varImp(rf_fit_neg1)
imp_neg1 <- imp_neg1$importance
imp_neg1 <- rownames_to_column(imp_neg1, "mz")
imp_neg1$mz <- as.numeric(gsub("X", "", imp_neg1$mz))
imp_neg1 <- imp_neg1[order(-imp_neg1$Overall),]

imp_neg2 <- varImp(rf_fit_neg2)
imp_neg2 <- imp_neg2$importance
imp_neg2 <- rownames_to_column(imp_neg2, "mz")
imp_neg2$mz <- as.numeric(gsub("X", "", imp_neg2$mz))
imp_neg2 <- imp_neg2[order(-imp_neg2$Overall),]

imp_neg3 <- varImp(rf_fit_neg3)
imp_neg3 <- imp_neg3$importance
imp_neg3 <- rownames_to_column(imp_neg3, "mz")
imp_neg3$mz <- as.numeric(gsub("X", "", imp_neg3$mz))
imp_neg3 <- imp_neg3[order(-imp_neg3$Overall),]

imp_neg4 <- varImp(rf_fit_neg4)
imp_neg4 <- imp_neg4$importance
imp_neg4 <- rownames_to_column(imp_neg4, "mz")
imp_neg4$mz <- as.numeric(gsub("X", "", imp_neg4$mz))
imp_neg4 <- imp_neg4[order(-imp_neg4$Overall),]

imp_neg5 <- varImp(rf_fit_neg5)
imp_neg5 <- imp_neg5$importance
imp_neg5 <- rownames_to_column(imp_neg5, "mz")
imp_neg5$mz <- as.numeric(gsub("X", "", imp_neg5$mz))
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
imp_pos1$mz <- as.numeric(gsub("X", "", imp_pos1$mz))
imp_pos1 <- imp_pos1[order(-imp_pos1$Overall),]

imp_pos2 <- varImp(rf_fit_pos2)
imp_pos2 <- imp_pos2$importance
imp_pos2 <- rownames_to_column(imp_pos2, "mz")
imp_pos2$mz <- as.numeric(gsub("X", "", imp_pos2$mz))
imp_pos2 <- imp_pos2[order(-imp_pos2$Overall),]

imp_pos3 <- varImp(rf_fit_pos3)
imp_pos3 <- imp_pos3$importance
imp_pos3 <- rownames_to_column(imp_pos3, "mz")
imp_pos3$mz <- as.numeric(gsub("X", "", imp_pos3$mz))
imp_pos3 <- imp_pos3[order(-imp_pos3$Overall),]

imp_pos4 <- varImp(rf_fit_pos4)
imp_pos4 <- imp_pos4$importance
imp_pos4 <- rownames_to_column(imp_pos4, "mz")
imp_pos4$mz <- as.numeric(gsub("X", "", imp_pos4$mz))
imp_pos4 <- imp_pos4[order(-imp_pos4$Overall),]

imp_pos5 <- varImp(rf_fit_pos5)
imp_pos5 <- imp_pos5$importance
imp_pos5 <- rownames_to_column(imp_pos5, "mz")
imp_pos5$mz <- as.numeric(gsub("X", "", imp_pos5$mz))
imp_pos5 <- imp_pos5[order(-imp_pos5$Overall),]

# isolate the top 100 most important from each model
imp_pos1 <- head(imp_pos1, n = 100)
imp_pos2 <- head(imp_pos2, n = 100)
imp_pos3 <- head(imp_pos3, n = 100)
imp_pos4 <- head(imp_pos4, n = 100)
imp_pos5 <- head(imp_pos5, n = 100)

imp_pos <- rbind(imp_pos1, imp_pos2, imp_pos3, imp_pos4, imp_pos5)

# set the phenotype to phenotype + a number so that it is unique and can be a row name
train_neg$phenotype <- paste(train_neg$phenotype, 1:66, sep = "_")
rownames(train_neg) <- NULL
train_neg <- column_to_rownames(train_neg, "phenotype")
# transpose and add an mz column
train_neg_t <- as.data.frame(t(train_neg))
train_neg_t <- rownames_to_column(train_neg_t, "mz")
train_neg_t$mz <- as.numeric(gsub("X", "", train_neg_t$mz))

imp_var_neg <- train_neg_t %>%
  filter_all(any_vars(.%in% imp_neg$mz))
imp_var_neg <- column_to_rownames(imp_var_neg, "mz")
imp_var_neg_t <- as.data.frame(t(imp_var_neg))
imp_var_neg_t <- rownames_to_column(imp_var_neg_t, "phenotype")
imp_var_neg_t$phenotype <- as.factor(ifelse(grepl(imp_var_neg_t$phenotype,
                                                 pattern = "differentiated"),
                                            "differentiated", "naive"))
train_neg <- imp_var_neg_t
# 201 peaks

# set the phenotype to phenotype + a number so that it is unique and can be a row name
train_pos$phenotype <- paste(train_pos$phenotype, 1:66, sep = "_")
rownames(train_pos) <- NULL
train_pos <- column_to_rownames(train_pos, "phenotype")
# transpose and add an mz column
train_pos_t <- as.data.frame(t(train_pos))
train_pos_t <- rownames_to_column(train_pos_t, "mz")
train_pos_t$mz <- as.numeric(gsub("X", "", train_pos_t$mz))

imp_var_pos <- train_pos_t %>%
  filter_all(any_vars(.%in% imp_pos$mz))
imp_var_pos <- column_to_rownames(imp_var_pos, "mz")
imp_var_pos_t <- as.data.frame(t(imp_var_pos))
imp_var_pos_t <- rownames_to_column(imp_var_pos_t, "phenotype")
imp_var_pos_t$phenotype <- as.factor(ifelse(grepl(imp_var_pos_t$phenotype,
                                                  pattern = "differentiated"),
                                            "differentiated", "naive"))
train_pos <- imp_var_pos_t
# 241 peaks

test_neg$phenotype <- paste(test_neg$phenotype, 1:16, sep = "_")
rownames(test_neg) <- NULL
test_neg <- column_to_rownames(test_neg, "phenotype")
# transpose and add an mz column
test_neg_t <- as.data.frame(t(test_neg))
test_neg_t <- rownames_to_column(test_neg_t, "mz")
test_neg_t$mz <- as.numeric(gsub("X", "", test_neg_t$mz))

test_neg_t <- test_neg_t %>%
  filter_all(any_vars(.%in% imp_neg$mz))
test_neg_t <- column_to_rownames(test_neg_t, "mz")
test_neg <- as.data.frame(t(test_neg_t))
test_neg <- rownames_to_column(test_neg, "phenotype")
test_neg$phenotype <- as.factor(ifelse(grepl(test_neg$phenotype,
                                             pattern = "differentiated"),
                                       "differentiated", "naive"))

colnames(train_neg) <- as.character(colnames(train_neg))

test_pos$phenotype <- paste(test_pos$phenotype, 1:16, sep = "_")
rownames(test_pos) <- NULL
test_pos <- column_to_rownames(test_pos, "phenotype")
# transpose and add an mz column
test_pos_t <- as.data.frame(t(test_pos))
test_pos_t <- rownames_to_column(test_pos_t, "mz")
test_pos_t$mz <- as.numeric(gsub("X", "", test_pos_t$mz))

test_pos_t <- test_pos_t %>%
  filter_all(any_vars(.%in% imp_pos$mz))
test_pos_t <- column_to_rownames(test_pos_t, "mz")
test_pos <- as.data.frame(t(test_pos_t))
test_pos <- rownames_to_column(test_pos, "phenotype")
test_pos$phenotype <- as.factor(ifelse(grepl(test_pos$phenotype,
                                             pattern = "differentiated"),
                                       "differentiated", "naive"))

colnames(train_pos) <- as.character(colnames(train_pos))

set.seed(2017)
model_neg <- randomForest(x = train_neg[-1],
                          y = train_neg$phenotype,
                          data = train_neg,
                          ntree = 500,
                          mtry = 20,
                          importance = TRUE,
                          proximity = TRUE)

print(model_neg)
plot(model_neg)

model_pos <- randomForest(x = train_pos[-1],
                          y = train_pos$phenotype,
                          data = train_pos,
                          ntree = 500,
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
for (i in c(1:200)){
  set.seed(1230)
  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
                      y = train_neg$phenotype,
                      data = train_neg,
                      ntree = 300,
                      mtry = i,
                      importance = TRUE)

  modellist_neg$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}

# mtry =3  has the lowest error rate, 201 peaks, error rate of 6%

modellist_pos <- data.frame(mtry = 1:240)
for (i in c(1:240)){
  set.seed(1230)
  fit <- randomForest(x = (train_neg %>% dplyr :: select(-phenotype)),
                      y = train_neg$phenotype,
                      data = train_neg,
                      ntree = 300,
                      mtry = i,
                      importance = TRUE)

  modellist_pos$err_rate[i] <- fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}
# mtry =   has the lowest error rate, 240 peaks, error rate of 6%


# permutation of data
permuted_neg <- rbind(train_neg, test_neg)
x <- c("differentiated", "naive")
# get all permutations
set.seed(1230)
false_samples_neg <- sample(x, 82, replace = TRUE)

permuted_neg$phenotype <- as.factor(false_samples_neg)
permuted_predict_neg <- predict(model_neg, newdata = permuted_neg[-1])
confusionMatrix(permuted_predict_neg, permuted_neg$phenotype)
# accuracy 0.5122
results_permuted_neg <- as.data.frame(cbind("actual" = permuted_neg$phenotype,
                                            "prediction" = permuted_predict_neg))
roc_permuted_neg <- roc(results_permuted_neg$actual, results_permuted_neg$prediction)
auc_permuted_neg <- auc(roc_permuted_neg)

# permutation of data
permuted_pos <- rbind(train_pos, test_pos)
x <- c("differentiated", "naive")
# get all permutations
set.seed(1230)
false_samples_pos <- sample(x, 82, replace = TRUE)

permuted_pos$phenotype <- as.factor(false_samples_pos)
permuted_predict_pos <- predict(model_pos, newdata = permuted_pos[-1])
confusionMatrix(permuted_predict_pos, permuted_pos$phenotype)
# accuracy 0.5122
results_permuted_pos <- as.data.frame(cbind("actual" = permuted_pos$phenotype,
                                            "prediction" = permuted_predict_pos))
roc_permuted_pos <- roc(results_permuted_pos$actual, results_permuted_pos$prediction)
auc_permuted_pos <- auc(roc_permuted_pos)

# AUC and ROC, first prep results from above models
test_neg1$phenotype <- as.factor(ifelse(grepl(test_neg1$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_neg2$phenotype <- as.factor(ifelse(grepl(test_neg2$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_neg3$phenotype <- as.factor(ifelse(grepl(test_neg3$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_neg4$phenotype <- as.factor(ifelse(grepl(test_neg4$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_neg5$phenotype <- as.factor(ifelse(grepl(test_neg5$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))

test_pos1$phenotype <- as.factor(ifelse(grepl(test_pos1$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_pos2$phenotype <- as.factor(ifelse(grepl(test_pos2$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_pos3$phenotype <- as.factor(ifelse(grepl(test_pos3$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_pos4$phenotype <- as.factor(ifelse(grepl(test_pos4$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))
test_pos5$phenotype <- as.factor(ifelse(grepl(test_pos5$phenotype,
                                              pattern = "differentiated"),
                                        "differentiated", "naive"))

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

# ROC curve of each fold/section
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
  scale_color_manual(breaks = c("differentiated", "naive"),
                     values = c("red4", "deepskyblue2")) +
  ggtitle("MDS plot (neg) using random forest proximities")

ggplot(mds_pos_data, aes(x = x, y = y, label = Sample)) +
  geom_point(aes(color = Phenotype), size = 3) +
  theme_classic() +
  theme(legend.position = "bottom",) +
  theme(text = element_text(size = 18)) +
  xlab(paste("MDS1 - ", mds_var_pos[1], "%", sep = "")) +
  ylab(paste("MDS2 - ", mds_var_pos[2], "%", sep = "")) +
  scale_color_manual(breaks = c("differentiated", "naive"),
                     values = c("red4", "deepskyblue2")) +
  ggtitle("MDS plot (pos) using random forest proximities")
