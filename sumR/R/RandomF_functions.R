## data should be in transposed form (mz as columns and samples column named =samples )
#' @title Redundant variables removal
#' @description This function removes correlated peaks at a certain cutoff value
#' @param data the dataframe with m/z as columns and a sample column named samples
#' @param cuttoff a certain value of cuttoff for high correlation (default=0.75)
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom caret findCorrelation
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom stats cor
redvar_removal <- function (data,cutoff=0.75){
  data<-column_to_rownames(data,"mz")
  data<-as.data.frame(t(data))
  cor_matrix <- cor(data)
  high_cor <- findCorrelation(cor_matrix, cutoff=cutoff)
  high_cor_removal<-data[,high_cor] ## highly correlated peaks
  data<- data %>% dplyr :: select(-colnames(high_cor_removal)) ## removing of high cor peaks
  data<-rownames_to_column(data,"samples")
  return(data)
}



## feature selection
#' @title feature selection
#' @description This function select the desired numbers of most imp features for random forest model
#' @param data_t the dataframe with m/z as columns and a sample column named samples
#' @param subsets how many peaks do u want to be selected or the range of peaks
#' @param seed global seed for reproducible results
#' @param method default is "repeatedcv" the method for feature selection is repeated cross validation
#' @param repeats numbers of repeats for cross validation default = 5
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @importFrom caret rfeControl rfe rfFuncs
#' @importFrom dplyr select
#' @importFrom magrittr %>%
imp_select<-function(data,data_t,subsets,seed,method="repeatedcv",repeats=5,class1,class2){
  data_t$samples <- as.factor(ifelse(grepl(data_t$samples, pattern = class1), class1, class2))
  control <- rfeControl(functions = rfFuncs,method =method,repeats =repeats)
  set.seed(seed)
  feature_select_cor <- rfe(data_t %>% dplyr::select(-samples),
                            data_t$samples, rfeControl=control,sizes =subsets)
  return("feature_select_model"=feature_select_cor)
}

#' @title feature selection dataframe
#' @param data the master dataframe
#' @param feature_select_model feature select model from imp_select function
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @importFrom dplyr any_vars filter_all
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom caret predictors
#' @importFrom magrittr %>%
imp_data<-function(feature_select_model,data,class1,class2){
  data<-data %>% filter_all(any_vars(. %in% predictors(feature_select_model)))
  data<-column_to_rownames(data,"mz")
  data<-as.data.frame(t(data))
  data<-rownames_to_column(data,"samples")
  data$samples <- as.factor(ifelse(grepl(data$samples, pattern = class1), class1, class2))
  colnames(data)<-as.character(colnames(data))
  return(data)
}

#' @title feature selection plot
#' @param feature_select_model feature select model from imp_select function
#' @importFrom ggplot2 ggplot
plot_imp_select<-function(feature_select_model){
  plot<-ggplot(data = feature_select_model, metric = "Accuracy") + theme_bw()
  return(plot)
}


#' @title datasplit
#'
#' @param data the dataframe with m/z as columns and a sample column named samples
#' @param split_ratio 
#' @param seed global seed
#'
#' @importFrom caTools sample.split
data_split<-function(seed,data,split_ratio=0.8){
  set.seed(seed) ## reproducible results 
  split<- sample.split(data$samples, SplitRatio = split_ratio)
  training_set<-subset(data, split == TRUE)
  test_set<-subset(data, split == FALSE)
  return(list("training_set"=training_set,"test_set"=test_set))
}
# random forest cross validation ------------------------------------------


## Random Forest subsets splitting
#' @title Random Forest subsets splitting
#' @description This function split the main training set into multiple subsets acoording to K numbers for cross validation
#' @param k number of subsets splitted
#' @param data the training set
#' @param seed global seed for reproducible results
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
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


# Cross validation

#' @title feature selection by cross validation
#' @description This function apply random forest model to the splitted subsets and extract the important variable from each
#' model and the roc and auc for each model 
#' @param training_set the training set
#' @param test_set the test set
#' @param k number of subsets splitted 
#' @param seed global seed for reproducible results
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @param n number of imp features to be selected in each model
#' @importFrom randomForest randomForest importance
#' @importFrom caret rfe confusionMatrix
#' @importFrom dplyr select arrange desc
#' @importFrom pROC roc auc
#' @importFrom magrittr %>%
#' @importFrom utils head
#' @importFrom stats predict
randomForest_CV<- function (training_set,test_set,class1,class2,k,seed,n){
  rflist <- rf_group(k = k,data=training_set,class1,class2, seed = seed)  
  data<-training_set
  cfmatrix <- list() # store the confusion matrix table of each group
  auc_all <- list()  # store the auc of each group
  roc_all <- list() # store the ROC curve of each group
  cfmatrix_summary <- list() # store the confusion matrix summary of each group
  imp<-list()
  imp_all<-data.frame()
  for (i in 1:k){
    train <- data[-rflist[[i]],]  
    test <- data[rflist[[i]],]
    set.seed(seed)
    model<-randomForest(x=(train %>% dplyr::select(-samples)),
                        y=train$samples,data=train,ntree=500,importance = T,proximity = T)
    tmpimp <- as.data.frame(importance(model))
    tmpimp<-rownames_to_column(tmpimp,"mz")
    tmpimp <- tmpimp %>% arrange(desc(MeanDecreaseAccuracy))
    imp [[i]] <-  tmpimp
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
    imp[i]<-head(imp[[i]],n=n)
    imp_all<-rbind(imp_all,imp[i])
  }
  colnames(imp_all)<-"mz"
  return(list("auc_all"=auc_all,"roc_all"=roc_all,"imp"=imp,
              "cfmatrix_summary"=cfmatrix_summary,"imp_all"=imp_all))
}


#' @title Cross-validation data
#' @param data the training set or test set
#' @param imp_all the dataframe of all imp m/z from the list of randomForest_CV function
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @importFrom dplyr filter_all any_vars
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom magrittr %>%
CV_data<-function(data,imp_all,class1,class2){
  data$samples<-paste(data$samples,1:nrow(data),sep="_") ## we just do this step so we can move the samples to row names .. 
  rownames(data)<-NULL
  data<-column_to_rownames(data,"samples")
  data_t<-as.data.frame(t(data))
  data_t<-rownames_to_column(data_t,"mz")
  imp_var<-data_t %>% filter_all(any_vars(. %in% imp_all$mz)) ## obtaining the intensities of the all imp mz 
  imp_var<-column_to_rownames(imp_var,"mz")
  imp_var_t<-as.data.frame(t(imp_var))
  imp_var_t<-rownames_to_column(imp_var_t,"samples")
  imp_var_t$samples<-as.factor(ifelse(grepl(imp_var_t$samples, pattern = class1), class1, class2))
  data<-imp_var_t
  colnames(data)<-as.character(colnames(data))
  return(data)
}

# the final model

#' @title random forest model
#' @description This function make the final random forest model with the selected features  
#' @param training_set the training set
#' @param test_set the test set
#' @param mtry number of mtry for the model
#' @param seed global seed for reproducible results
#' @param ntree number of trees for the model 
#' @importFrom randomForest randomForest
#' @importFrom dplyr select
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom magrittr %>%
#' @importFrom stats predict
RF_model<-function(training_set,test_set,mtry,ntree,seed){
  ## the final model
  set.seed(seed) 
  model<-randomForest(x=(training_set %>% dplyr::select(-samples)),y=training_set$samples,data=training_set,ntree=ntree,mtry=mtry,
                      importance = T,proximity = T) #importance=T so we can plot imp var 
  test_set$samples <- as.factor(test_set$samples)
  prediction <- predict(model,newdata=test_set,type="class")
  confusion_matrix<-confusionMatrix(table(data=prediction, reference = test_set$samples))
  results<-as.data.frame(cbind("Actual"=test_set$samples,"Prediction"=prediction))
  rocfinal<- roc(results$Actual,results$Prediction )
  aucfinal<- auc(rocfinal)
  return(list("results"=results,"rocfinal"=rocfinal,"aucfinal"=aucfinal,"model"=model,"confusion_matrix"=confusion_matrix))
}


## choosing the best mtry


#' @title choose mtry
#' @description This function choose the mtry according to the user preference (Accuracy, Sensitivity, Specificity)
#' @param training_set the training set
#' @param test_set the test set
#' @param seed global seed for reproducible results
#' @param ntree number of trees used in the original model
#' @importFrom randomForest randomForest
#' @importFrom dplyr select
#' @importFrom caret confusionMatrix
#' @importFrom magrittr %>%
#' @importFrom stats predict
mtry_select<-function(training_set,test_set,seed,ntree){
  model_list <- data.frame(mtry=1:(length(training_set)-1))
  for (i in c(1:(length(training_set)-1))) {
    set.seed(seed)
    fit <- randomForest(x=(training_set %>% dplyr::select(-samples)),
                        y=training_set$samples,data=training_set,ntree=ntree,mtry=i,importance = T)
    model_list$err_rate[i]<-fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
    ## prediction of model
    test_set$samples <- as.factor(test_set$samples)
    prediction <- predict(fit,newdata=test_set,type="class")
    cf<-confusionMatrix(data=prediction, reference = test_set$samples)
    model_list$Accuracy[i]<-cf[["overall"]][["Accuracy"]]
    model_list$Sensitivity[i]<-cf[["byClass"]][["Sensitivity"]]
    model_list$Specificity[i]<-cf[["byClass"]][["Specificity"]]
  }
  return(model_list)
}




# permutation of data  ----------------------------------------------------



#' @title data permutation
#' @description This function apply data permutation to test the accuracy of the model against permutated data  
#' @param training_set the training set
#' @param test_set the test set
#' @param model the random forest model  
#' @param seed global seed for reproducible results
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @param n number permutation
#' @param sample_no the number of all samples(training+test set)
#' @importFrom  dplyr select
#' @importFrom caret confusionMatrix
#' @importFrom magrittr %>%
#' @importFrom graphics hist
#' @importFrom stats predict
permutation.test <- function(training_set,test_set, model, n,class1,class2,sample_no){
  permuted_data<-rbind(training_set,test_set)
  permuted_df <- data.frame(n=1:n)  
  for (i in c(1:n)) {
    x <- c(class1,class2)
    false_samples<-sample(x, sample_no, replace = TRUE)
    permuted_data$samples<-as.factor(false_samples)
    permuted_predict<-predict(model,newdata = permuted_data %>% dplyr::select(-samples))
    cf<-confusionMatrix(permuted_predict,permuted_data$samples)
    cf<-as.data.frame(cf$overall)
    permuted_df$accuracy[i]<-cf[1,]
  }
  plot<-hist(permuted_df$accuracy)
  print(plot)
  return(list(plot,permuted_df))
}
