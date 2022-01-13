#
#' @title Redundant variables removal
#' @description This function removes correlated peaks at a certain cutoff value
#' @param data the dataframe with mz column and samples columns
#' @param cuttoff a certain value of cuttoff for high correlation
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom caret findCorrelation
#' @importFrom dplyr select
redvar_removal <- function (data,cutoff){
data <- column_to_rownames(data,"mz")
  data <- as.data.frame(t(data))

  cor_matrix <- cor(data)

  # find attributes that are highly corrected (ideally >0.75)
  high_cor <- findCorrelation(cor_matrix, cutoff=cutoff)

  high_cor_removal<-data[,high_cor] ## highly correlated peaks
  data<- data %>% dplyr::select(-colnames(high_cor_removal)) ## removing of high cor peaks
  data<-rownames_to_column(data,"samples")

  return(data)
}





## feature selection
#' @title feature selection
#' @description This function select the desired numbers of most imp features for random forest model
#' @param data_t the dataframe with m/z as columns and a sample column named samples
#' @param data the master dataframe
#' @param subsets how many peaks do u want to be selected or the range of peaks
#' @param seed global seed for reproducible results
#' @param method "repeatedcv" the method for feature selection is repeated cross validation
#' @param repeats numbers of repeats for cross validation
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom caret rfeControl
#' @importFrom caret rfe
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr filter_all
imp_select<-function(data,data_t,subsets,seed,method,repeats,class1,class2){


  data_t$samples <- as.factor(ifelse(grepl(data_t$samples, pattern = class1), class1, class2))

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
  colnames(data)<-as.character(colnames(data))

  return(data)
}


#' @title datasplit
#' @param data the dataframe with m/z as columns and a sample column named samples
#' @param split ratio a certain value for split ratio
#' @param seed global seed
#' @importFrom caTools sample.split
data_split<-function(seed,data,split_ratio){
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

#' @title feature selection
#' @description This function apply random forest model to the splitted subsets and extract the important variable from each
#' model and the roc and auc for each model
#' @param training_set the training set
#' @param test_set the test set
#' @param k number of subsets splitted
#' @param seed global seed for reproducible results
#' @param class1 sample group 1 name for renaming
#' @param class2 sample group 2 name for renaming
#' @param n number of imp features to be selected in each model
#' @importFrom randomForest randomForest
#' @importFrom randomForest importance
#' @importFrom dplyr select
#' @importFrom caret rfe
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr filter_all
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc
#' @importFrom pROC auc
randomForest_CV<- function (training_set,test_set,class1,class2,k,seed,n){



  rflist <- rf_group(k = k,data=training_set,class1,class2, seed = seed)
  data<-training_set

  cfmatrix <- list() # store the confusion matrix table of each group
  auc_all <- list()  # store the auc of each group
  roc_all <- list() # store the ROC curve of each group
  cfmatrix_summary <- list() # store the confusion matrix summary of each group
  imp<-list()
  set.seed(seed)
  for (i in 1:k){
    train <- data[-rflist[[i]],]
    test <- data[rflist[[i]],]
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



  }


  imp_1<-head(imp[[1]],n=n)
  imp_2<-head(imp[[2]],n=n)
  imp_3<-head(imp[[3]],n=n)
  imp_4<-head(imp[[4]],n=n)
  imp_5<-head(imp[[5]],n=n)

  imp_all<-rbind(imp_1,imp_2,imp_3,imp_4,imp_5)

  training_set$samples<-paste(training_set$samples,1:nrow(training_set),sep="_") ## we just do this step so we can move the samples to row names ..
  rownames(training_set)<-NULL
  training_set<-column_to_rownames(training_set,"samples")

  training_set_t<-as.data.frame(t(training_set))

  training_set_t<-rownames_to_column(training_set_t,"mz")

  imp_var<-training_set_t %>% filter_all(any_vars(. %in% imp_all$mz)) ## obtaining the intensities of the all imp mz
  imp_var<-column_to_rownames(imp_var,"mz")
  imp_var_t<-as.data.frame(t(imp_var))
  imp_var_t<-rownames_to_column(imp_var_t,"samples")
  imp_var_t$samples<-as.factor(ifelse(grepl(imp_var_t$samples, pattern = class1), class1, class2))
  training_set<-imp_var_t


  ## for test set ## what we do here is we make the test set with same selected imp peaks so we can test the final model


  test_set$samples<-paste(test_set$samples,1:nrow(test_set),sep="_")
  rownames(test_set)<-NULL
  test_set<-column_to_rownames(test_set,"samples")

  test_set_t<-as.data.frame(t(test_set))

  test_set_t<-rownames_to_column(test_set_t,"mz")

  test_set_t<-test_set_t %>% filter_all(any_vars(. %in% imp_all$mz))
  test_set_t<-column_to_rownames(test_set_t,"mz")
  test_set<-as.data.frame(t(test_set_t))
  test_set<-rownames_to_column(test_set,"samples")
  test_set$samples<-as.factor(ifelse(grepl(test_set$samples, pattern = class1), class1, class2))



  colnames(training_set)<-as.character(colnames(training_set))

  return(list("training_set"=training_set,"test_set"=test_set,"auc_all"=auc_all,"roc_all"=roc_all,"imp"=imp,
              "cfmatrix_summary"=cfmatrix_summary))
}



## the final model

#' @title random forest model
#' @description This function make the final random forest model with the selected features
#' @param training_set the training set
#' @param test_set the test set
#' @param mtry number of mtry for the model
#' @param seed global seed for reproducible results
#' @param ntree number of trees for the model
#' @importFrom randomForest randomForest
#' @importFrom dplyr select
#' @importFrom caret rfe
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc
#' @importFrom pROC auc
RF_model<-function(training_set,test_set,mtry,ntree,seed){

  ## the final model
  set.seed(seed)
  model<-randomForest(x=(training_set %>% dplyr::select(-samples)),y=training_set$samples,data=training_set,ntree=ntree,mtry=mtry,
                      importance = T,proximity = T) #importance=T so we can plot imp var
  # proximity = T so we can plot proximity plot

  print(model)
  plot(model)


  ## prediction of model
  test_set$samples <- as.factor(test_set$samples)

  prediction <- predict(model,newdata=test_set,type="class")

  print(confusionMatrix(table(data=prediction, reference = test_set$samples)))

  results<-as.data.frame(cbind("Actual"=test_set$samples,"Prediction"=predtrain))

  rocfinal<- roc(results$Actual,results$Prediction )
  aucfinal<- auc(rocfinal)
  return(list("results"=results,"rocfinal"=rocfinal,"aucfinal"=aucfinal))
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
#' @importFrom lattice histogram
permutation.test <- function(training_set,test_set, model, n,class1,class2,sample_no){
  data<-rbind(training_set,test_set)

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
#get all permutations


