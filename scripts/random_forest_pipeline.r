library(plyr)
library(tidyverse)
library(tibble)
library(caTools) #forsplitting the dataset
library(randomForest) # random forest
library(caret) #feature selection 
library(pROC) ## roc curve 
library(rfUtilities) ## cross validation 

#  feature selection for random forest -----------------------------------------------------------

# master_random<-master_df_t
# 
# 
# ## removal of highly correlated peaks
# ### removing correlated data (redundant features)
# cor_matrixeat <- cor(master_random)
# 
# # find attributes that are highly corrected (ideally >0.75)
# high_cor_feat <- findCorrelation(cor_matrix_feat, cutoff=0.75)
# # print indexes of highly correlated attributes
# print(high_cor_feat)
# 
# high_cor_removal<-master_random[,high_cor_feat] ## highly correlated peaks
# master_random<- master_random %>% dplyr :: select(-colnames(high_cor_removal)) ## removing of high cor peaks
# 
# master_random<-rownames_to_column(master_random,"samples")
# 
## data should be in transposed form (mz as columns and samples column named =samples )
redvar_removal <- function (data,cutoff){
  data<-column_to_rownames(data,"samples")
  
  cor_matrix <- cor(data)
  
  # find attributes that are highly corrected (ideally >0.75)
  high_cor <- findCorrelation(cor_matrix, cutoff=cutoff)
  
  high_cor_removal<-data[,high_cor] ## highly correlated peaks
  data<- data %>% dplyr :: select(-colnames(high_cor_removal)) ## removing of high cor peaks
  data<-rownames_to_column(data,"samples")
  
  return(data)
}

test<-redvar_removal(sig_data_t,0.75)

## renaming the samples so they all have the same class name



# subsets <- c(150:300,500) ##choosing a subsets we need for feature selections .. we can either put the whole range of peaks
# ## or pre-set a certain number to be selected from like 50,100,150 .. etc or (100:500) or we can mix both
# ## its a rule of thumb that the feature has to be max three times the sample number
# ## so if we have like 50 sample our peaks should be around 150
# 
# ## control for our feature selection .. is using a random forest algorithm with repeated cross validation 
# control <- rfeControl(functions = rfFuncs,method ="repeatedcv",repeats =5)
# # run the RFE algorithm
# set.seed(1e6)
# 
# feature_select <- rfe(master_random %>% dplyr::select(-samples),
#                           master_random$samples, rfeControl=control,sizes =subsets)
# ## now we run our feature selection
# # summarize the results
# print(feature_select)
# 
# # we can see here the number of selected peaks
# # plot the results
# ggplot(data = feature_select, metric = "Accuracy") + theme_bw() 
# ggplot(data = feature_select, metric = "Kappa") + theme_bw()
# 
# mz_select<-master_df_pos
# 
# ## extracting our selected features and obtaining their scaled intensities to be used in the random forest 
# mz_select<-mz_select %>% filter_all(any_vars(. %in% predictors(feature_select)))
# 
# mz_select<-column_to_rownames(mz_select,"mz")
# mz_select_t<-as.data.frame(t(mz_select))
# mz_select_t<-rownames_to_column(mz_select_t,"samples")

imp_select<-function(data,data_t,subsets,seed,method,repeats,class1,class2){
  
  
  data_t$samples <- as.factor(ifelse(grepl(data_t$samples, pattern = "X"), "X", "Y"))
  
  control <- rfeControl(functions = rfFuncs,method =method,repeats =repeats)
  # run the RFE algorithm
  set.seed(seed)
  
  feature_select_cor <- rfe(data_t %>% dplyr::select(-samples),
                            data_t$samples, rfeControl=control,sizes =subsets)
  
  # summarize the results
  print(feature_select_cor)
  
  # plot the results
  
  plot<-ggplot(data = feature_select_cor, metric = "Accuracy") + theme_bw() 
  print(plot)
  
  
  data<-data %>% filter_all(any_vars(. %in% predictors(feature_select_cor)))
  
  data<-column_to_rownames(data,"mz")
  data<-as.data.frame(t(data))
  data<-rownames_to_column(data,"samples")
  data$samples <- as.factor(ifelse(grepl(data$samples, pattern = class1), class1, class2))
  
  return(data)
}

mz_select_t<-imp_select(master_df_pos,sig_data_t,50,1e6,"repeatedcv",2,"X","Y")


colnames(mz_select_t)<-as.character(colnames(mz_select_t))


# random forest cross validation ------------------------------------------

# The concept is to split your data set into primary training set and primary test set
# Then split your primary training set into secondary training set and secondary test set and use them to make your first random forest model 
# Repeat this step 5 times 
# Extract the import variable from each model
# Use these important variable to build the final model and test this model using the primary test set we set aside from the beginning 

set.seed(1e6) ## reproducible results 
split<- sample.split(mz_select_t$samples, SplitRatio = 0.8)

training_set<-subset(mz_select_t, split == TRUE)
test_set<-subset(mz_select_t, split == FALSE)

# # splitting again for testing - first test
# set.seed(1e6) 
# split_f<- sample.split(training_set$samples, SplitRatio = 0.8)
# 
# training_set_f<-subset(training_set_f, split_f_1 == TRUE)
# test_set_f<-subset(training_set_f, split_f_1 == FALSE)
# 
# colnames(training_set_f)<-as.character(colnames(training_set_f))
# 
# set.seed(1e6) ## the first model
# model_1<-randomForest(x=training_set_1[-1],y=training_set_1$samples,data=training_set,ntree=500)
# print(model_1)
# plot(model_1)
# predtrain_1 <- predict(model_1,newdata=test_set_1)
# 
# 
# actual_sample_1 <- as.factor(test_set_1$samples)
# confusionMatrix(data=predtrain_1, reference = actual_sample_1)
# results_1<-as.data.frame(cbind("actual"=test_set_1$samples,"prediction"=predtrain_1))
# 
# 
# imp_1<-varImp(model_1)
# imp_1<-rownames_to_column(imp_1,"mz")
# imp_1<-imp_1[order(-imp_1$Overall),] ## order the imp variables in descending order 
# 
# # splitting again for testing - second test .. repeating as the first model and so on
# set.seed(1e6) 
# split_2<- sample.split(training_set$samples, SplitRatio = 0.8)
# training_set_2<-subset(training_set, split_2 == TRUE)
# test_set_2<-subset(training_set, split_2 == FALSE)
# 
# colnames(training_set_2)<-as.character(colnames(training_set_2))
# set.seed(1e6) 
# model_2<-randomForest(x=training_set_2[-1],y=training_set_2$samples,data=training_set,ntree=500)
# print(model_2)
# plot(model_2)
# predtrain_2 <- predict(model_2,newdata=test_set_2)
# 
# 
# actual_sample_2 <- as.factor(test_set_2$samples)
# confusionMatrix(data=predtrain_2, reference = actual_sample_2)
# results_2<-as.data.frame(cbind("actual"=test_set_2$samples,"prediction"=predtrain_2))
# 
# 
# imp_2<-varImp(model_2)
# imp_2<-rownames_to_column(imp_2,"mz")
# imp_2<-imp_2[order(-imp_2$Overall),]
# 
# 
# # splitting again for testing - third test
# set.seed(1e6) 
# split_3<- sample.split(training_set$samples, SplitRatio = 0.8)
# 
# training_set_3<-subset(training_set, split_3 == TRUE)
# test_set_3<-subset(training_set, split_3 == FALSE)
# 
# colnames(training_set_3)<-as.character(colnames(training_set_3))
# set.seed(1e6) 
# model_3<-randomForest(x=training_set_3[-1],y=training_set_3$samples,data=training_set,ntree=500)
# print(model_3)
# plot(model_3)
# predtrain_3 <- predict(model_3,newdata=test_set_3)
# 
# 
# actual_sample_3 <- as.factor(test_set_3$samples)
# confusionMatrix(data=predtrain_3, reference = actual_sample_3)
# results_3<-as.data.frame(cbind("actual"=test_set_3$samples,"prediction"=predtrain_3))
# 
# imp_3<-varImp(model_3)
# imp_3<-rownames_to_column(imp_3,"mz")
# imp_3<-imp_3[order(-imp_3$Overall),]
# 
# # splitting again for testing - fourth test
# set.seed(1e6) 
# split_4<- sample.split(training_set$samples, SplitRatio = 0.8)
# 
# training_set_4<-subset(training_set, split_4 == TRUE)
# test_set_4<-subset(training_set, split_4 == FALSE)
# 
# colnames(training_set_4)<-as.character(colnames(training_set_4))
# set.seed(1e6) 
# model_4<-randomForest(x=training_set_4[-1],y=training_set_4$samples,data=training_set,ntree=500)
# print(model_4)
# plot(model_4)
# predtrain_4 <- predict(model_4,newdata=test_set_4)
# 
# 
# actual_sample_4 <- as.factor(test_set_4$samples)
# confusionMatrix(data=predtrain_4, reference = actual_sample_4)
# results_4<-as.data.frame(cbind("actual"=test_set_4$samples,"prediction"=predtrain_4))
# 
# imp_4<-varImp(model_4)
# imp_4<-rownames_to_column(imp_4,"mz")
# imp_4<-imp_4[order(-imp_4$Overall),]
# 
# # splitting again for testing - fifth test
# set.seed(1e6) 
# split_5<- sample.split(training_set$samples, SplitRatio = 0.8)
# 
# training_set_5<-subset(training_set, split_5 == TRUE)
# test_set_5<-subset(training_set, split_5 == FALSE)
# colnames(training_set_5)<-as.character(colnames(training_set_5))
# set.seed(1e6) 
# 
# model_5<-randomForest(x=training_set_5[-1],y=training_set_5$samples,data=training_set,ntree=500)
# print(model_5)
# plot(model_5)
# predtrain_5 <- predict(model_5,newdata=test_set_5)
# 
# 
# actual_sample_5 <- as.factor(test_set_5$samples)
# confusionMatrix(data=predtrain_5, reference = actual_sample_5)
# results_5<-as.data.frame(cbind("actual"=test_set_5$samples,"prediction"=predtrain_5))
# 
# imp_5<-varImp(model_5)
# imp_5<-rownames_to_column(imp_5,"mz")
# imp_5<-imp_5[order(-imp_5$Overall),]
# 
# ### the top 100 imp var in each model
# imp_1<-head(imp_1,n=100)
# imp_2<-head(imp_2,n=100)
# imp_3<-head(imp_3,n=100)
# imp_4<-head(imp_4,n=100)
# imp_5<-head(imp_5,n=100)
# 
# imp<-rbind(imp_1,imp_2,imp_3,imp_4,imp_5)

# Cross validation
# assign samples to different group
rf_group <- function(k, data, class1 , class2 , seed){ 
  # k - how many parts will you divide your data into
  rflist <- list()
  set.seed(seed)
  
  # to keep a fixed ratio for X & Y in each group
  n_class1 <- sum(grepl(data$samples, pattern = class1))
  n_class2 <- sum(grepl(data$samples, pattern = class2))
  
  n1 <- rep(1:k,ceiling(n_class1/k))[1:n_class1]  
  n2 <- rep(1:k,ceiling(n_class2/k))[1:n_class2] 
  
  temp1 <- sample(n1,n_class1) # randomize n1
  temp2 <- sample(n2,n_class2) # randomize n2
  
  x <- 1:k
  dataseq1 <- 1:n_class1
  rflist1 <- lapply(x,function(x) dataseq1[temp1==x])  
  
  dataseq2 <- (n_class1+1):(n_class1+n_class2)
  rflist2 <- lapply(x,function(x) dataseq2[temp2==x])  
  
  for (i in 1:k){
    rflist[[i]] <- c(rflist1[[i]],rflist2[[i]])
  }
  
  return(rflist)
}

k <- 5
rflist <- list()

# each element of the list is a cross validation group
rflist <- rf_group(k = k,data=training_set,"X","Y", seed = 1e6)
rflist

data <- training_set
imp <- data.frame(MDG = rep(0, ncol(data)-1))   #store the prediction results

cfmatrix <- list() # store the confusion matrix table of each group
auc_all <- list()  # store the auc of each group
roc_all <- list() # store the ROC curve of each group
cfmatrix_summary <- list() # store the confusion matrix summary of each group

set.seed(1e6)
for (i in 1:k){
  train <- data[-rflist[[i]],]  
  test <- data[rflist[[i]],]
  model<-randomForest(x=(train %>% dplyr::select(-samples)),
                      y=train$samples,data=train,ntree=500,importance = T,proximity = T)
  
  #model <-randomForest(train[, -ncol(train)], train[, ncol(train)], proximity=TRUE,ntree = 500,importance = T)  
  tmpimp <- as.data.frame(importance(model))
  imp <- cbind(imp, tmpimp)
  mypred <- predict(model, newdata = test, type = "class")
  mypred<-as.data.frame(cbind("observed"=(test$samples),"predict"=mypred))
  
  xtab <- table(as.factor(mypred$predict), mypred$observed)
  print(confusionMatrix(xtab))
  cfmatrix_summary[[i]] <- confusionMatrix(xtab)
  cfmatrix[[i]] <- as.matrix(table(mypred$predict,mypred$observed))
  
  
  
  roc<- roc(mypred$observed,mypred$predict)
  auc<- auc(roc)
  
  auc_all[[i]] <- auc
  roc_all[[i]] <- roc
  
  
  
}



imp<-imp[,-1]
imp$mean <- apply(imp, 1, mean)
imp <- rownames_to_column(imp, "mz")
imp <- imp %>% arrange(desc(mean))

# collect the top 170 significant peaks information for final model annotation
imp <- imp[1:170, c(1,ncol(imp))]

training_set$samples<-paste(training_set$samples,1:40,sep="_") ## we just do this step so we can move the samples to row names .. 
#35 according to the samples number of training set
rownames(training_set)<-NULL
training_set<-column_to_rownames(training_set,"samples")

training_set_t<-as.data.frame(t(training_set))

training_set_t<-rownames_to_column(training_set_t,"mz")

imp_var<-training_set_t %>% filter_all(any_vars(. %in% imp$mz)) ## obtaining the intensities of the all imp mz 
imp_var<-column_to_rownames(imp_var,"mz")
imp_var_t<-as.data.frame(t(imp_var))
imp_var_t<-rownames_to_column(imp_var_t,"samples")
imp_var_t$samples<-as.factor(mgsub::mgsub(imp_var_t$samples,pattern,replacement))
training_set<-imp_var_t

## total 171 peaks after duplicates removal

## for test set ## what we do here is we make the test set with same selected imp peaks so we can test the final model 


test_set$samples<-paste(test_set$samples,1:10,sep="_")
rownames(test_set)<-NULL
test_set<-column_to_rownames(test_set,"samples")

test_set_t<-as.data.frame(t(test_set))

test_set_t<-rownames_to_column(test_set_t,"mz")

test_set_t<-test_set_t %>% filter_all(any_vars(. %in% imp$mz))
test_set_t<-column_to_rownames(test_set_t,"mz")
test_set<-as.data.frame(t(test_set_t))
test_set<-rownames_to_column(test_set,"samples")
test_set$samples<-as.factor(mgsub::mgsub(test_set$samples,pattern,replacement))

colnames(training_set)<-as.character(colnames(training_set))

## the final model
set.seed(1e6) 
model<-randomForest(x=training_set[-1],y=training_set$samples,data=training_set,ntree=500,mtry=23,
                      importance = T,proximity = T) #importance=T so we can plot imp var 
# proximity = T so we can plot proximity plot 

print(model)
plot(model)


## confusion matrix of model

pred_trainingseteat <- predict(model,newdata=training_set)

confusionMatrix(data=pred_trainingseteat, reference = training_set$samples)


## prediction of model 
test_set$samples <- as.factor(test_set$samples)

predtrain <- predict(model,newdata=test_set,type="class")

confusionMatrix(table(data=predtrain, reference = test_set$samples))



results<-as.data.frame(cbind("actual"=test_set$samples,"prediction"=predtrain))

## choosing the best mtry that we produce lowest err rate of the model 
modellist <- data.frame(mtry=1:171) ## 151 according to number of peaks
for (i in c(1:171)) {
  set.seed(1e6)
  fit <- randomForest(x=(training_set %>% dplyr::select(-samples)),
                      y=training_set$samples,data=training_set,ntree=500,mtry=i,importance = T)
  
  modellist$err_rate[i]<-fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}
## mtry=23 has the lowest error rate




# permutation of data  ----------------------------------------------------

## the concept here is we randomly change the samples name and see how the model will do when trying to classify the samples 
## if the accuracy is bad then our model is good 


permuted_data<-rbind(training_set,test_set)

permutation.test <- function(data, model, n,class1,class2,sample_no){
  permuted_df <- data.frame(n=1:n) 
  
  for (i in c(1:n)) {
    permuted_data<-data
    x <- c(class1,class2)
    false_samples<-sample(x, sample_no, replace = TRUE)
    permuted_data$samples<-as.factor(false_samples)
    permuted_predict<-predict(model,newdata = permuted_data %>% dplyr::select(-samples))
    cf<-confusionMatrix(permuted_predict,permuted_data$samples)
    cf<-as.data.frame(cf$overall)
    
    permuted_df$accuracy[i]<-cf[1,]
  }
  
  plot<-histogram(permuted_df$accuracy)
  print(plot)
  return(list(plot,permuted_df))
}
test1<-permutation.test(permuted_data, model_f, 1000,"X","Y",50)
#get all permutations


# ## Accuracy : 0.44
# results_permuted<-as.data.frame(cbind("actual"=permuted_data$samples,"prediction"=permuted_predict))
# roc_permuted<- roc(results_permuted$actual,results_permuted$prediction )
# auc_permuted<- auc(roc_permuted)

# AUC and roc ------------------------------------------------------------------





rocinal<- roc(results$actual,results$prediction )
aucinal<- auc(rocinal)

# ROC curve of each fold
par(mfrow=c(2,3))
plot(roc_all[[1]],col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_all[[1]], 3))))
plot(roc_all[[2]],col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_all[[2]], 3))))
plot(roc_all[[3]],col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_all[[3]], 3))))
plot(roc_all[[4]],col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_all[[4]], 3))))
plot(roc_all[[5]],col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_all[[5]], 3))))
plot(rocfinal,col = "blue", main = paste("final model, AUC:", as.character(round(aucfinal, 3))))



# Plot the important variables
varImpPlot(model,col="blue",pch= 2)


# MDS plot
distance_matrix <- as.dist(1-model$proximity)

mds <- cmdscale(distance_matrix, eig=TRUE, x.ret=TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
mds_var <- round(mds$eig/sum(mds$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds_values <- mds$points
mds_data <- data.frame(Sample=rownames(mds_values),
                       X=mds_values[,1],
                       Y=mds_values[,2],
                       Status=training_set$samples)

ggplot(data=mds_data, aes(x=X, y=Y, label=Sample)) + 
  geom_text(aes(color=Status)) +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("MDS plot using (1 - Random Forest Proximities)")



# model evaluation --------------------------------------------------------



rf_cv <- rf.crossValidation(model, training_set,seed=1e6, p=0.35, n=99, ntree=500) #p is the proportion of the data will be held for cross validation
print(rf_cv)
# Plot cross validation versus model producers accuracy
par(mfrow=c(1,2)) 
plot(rf_cv, type = "cv", main = "CV producers accuracy")
plot(rf_cv, type = "model", main = "Model producers accuracy")

# Plot cross validation versus model oob
par(mfrow=c(1,2)) 
plot(rf_cv, type = "cv", stat = "oob", main = "CV oob error")
plot(rf_cv, type = "model", stat = "oob", main = "Model oob error")
