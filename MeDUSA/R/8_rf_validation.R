# *** RandomForest Validate-----------------------------------------------------
#' rf_obj: Validate randomforest model
#'
#' @description
#' Random Forest is a robust tool for identifying the features that contribute
#' to correct phenotype prediction. For optimal performance of the model, it is
#' recommended to remove the highly correlated features before training and
#' testing the model. The mzlog_rf_correlation function utilizes the
#' findCorrelation function of the Caret package to isolate and remove the
#' highly correlated features. The new data set with these feature removed will
#' be used for the remaining functions of random forest. The mzlog_rf_select
#' function utilizes the Caret rfeControl function to determine which features
#' of the data set are the best predictors for phenotype and how many features
#' should be considered in the model. Using the reduced, likely predictors only
#' data set will reduce the resources needed for the remainder of random forest
#' processing. The mzlog_rf function divides the data set randomly into sections
#' for training and testing the model. The model will output a coefficient of
#' relevance for each feature, essentially ranking the features based on how
#' fundamental the feature was in predicting the phenotype. The rf_validate
#' function is used to compare the predictions of the model to the actual
#' phenotype of the sample. The accuracy of the model is also output.
#'
#' @param rf_obj \cr
#'   List : from mzlog_rf: list(model, test, train):  Expects "phenotype"
#' @param mtry_seed \cr
#'   List: which seed to use.
#' @param cores \cr
#'   Int: can I haz multithreading
#'
#'
#' @export
rf_validate <- function(rf_obj, mtry_range = c(1:200), trees = 500, mtry_seed = 1984, cores = 2){
  pred_train <- rf.predict(rf_obj$model, rf_obj$train)
  if( sum(!(pred_train$pred == rf_obj$train$phenotype)) ){
    warning("MeDUSA:: rf_validate: Train prediction did not match train model")
  }
  pred_test <- rf.predict(rf_obj$model, rf_obj$test)
  results <- as.data.frame(cbind("actual" = rf_obj$test$phenotype,
                                 "prediction" = pred_test$pred))

  ml <- data.frame(mtry = mtry_range)
  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    ml$err_rate <- pbapply::pbapply(ml, 1, rf.mtry_fit, cl = cl,
                                   data = rf_obj$train, seed = mtry_seed, trees = trees)
    best_mtry <- ml[ml$err_rate == min(ml$err_rate),]
    print(paste0("INFO:MeDUSA::rf_validate: Lowest error:", best_mtry$err_rate,
                " At mtry:", best_mtry$mtry ))
  }, finally={
   local.kill_threads(cl)
  })

  list(train_confusionMatrix <- pred_train$cm, test_confusionMatrix <- pred_test$cm, mtry_error <- ml)
}


# *** RandomForest Permuted -----------------------------------------------------
#' Rf_obj: Randomize phenotypes to negatively validate rf-model
#'
#' @description
#' Random Forest is a robust tool for identifying the features that contribute
#' to correct phenotype prediction. For optimal performance of the model, it is
#' recommended to remove the highly correlated features before training and
#' testing the model. The mzlog_rf_correlation function utilizes the
#' findCorrelation function of the Caret package to isolate and remove the
#' highly correlated features. The new data set with these feature removed will
#' be used for the remaining functions of random forest. The mzlog_rf_select
#' function utilizes the Caret rfeControl function to determine which features
#' of the data set are the best predictors for phenotype and how many features
#' should be considered in the model. Using the reduced, likely predictors only
#' data set will reduce the resources needed for the remainder of random forest
#' processing. The mzlog_rf function divides the data set randomly into sections
#' for training and testing the model. The model will output a coefficient of
#' relevance for each feature, essentially ranking the features based on how
#' fundamental the feature was in predicting the phenotype. The rf_validate
#' function is used to compare the predictions of the model to the actual
#' phenotype of the sample. The accuracy of the model is also output. To ensure
#' the accuracy of the model is truly based on biological phenotype, a
#' permutation of the data is used. The rf_permuted function randomly assigns a
#' phenotype to each of the samples and feeds this permuted data set through the
#' model. If the accuracy is significantly greater than 50%, the user should
#' reconsider the model.
#'
#' @param rf_obj \cr
#'   List : from mzlog_rf: list(model, test, train): Expects "phenotype"
#' @param mtry_seed \cr
#'   List: which seed to use.
#' @param cores \cr
#'   Int: can I haz multithreading
#'
#'
#' @export
rf_permuted <- function(rf_obj, pol, seed = 2540.632, plot = T){
  permuted <- rbind(rf_obj$train, rf_obj$test)
  attr <- as.vector(unique(permuted$phenotype))
  set.seed(seed)
  false_samples <- sample(attr, nrow(permuted), replace = TRUE)

  permuted$phenotype <- as.factor(false_samples)
  pred_permuted <- rf.predict(rf_obj$model, permuted)

  results_permuted<- as.data.frame(cbind("actual" = permuted$phenotype,
                                              "prediction" = pred_permuted$pred))
  roc_permuted <- pROC::roc(results_permuted$actual, results_permuted$prediction)
  auc_permuted <- pROC::auc(roc_permuted)

  plot_file <- paste0(local.output_dir(),local.dir_sep(),"RF_ROC_Permuted_",pol,".png")
  plot(roc_permuted, col = "Red", main = paste0(pol,"_RF_ROC Permuted"),
       sub = paste0("Acc:",pred_permuted$cm$overall["Accuracy"]," AUC:", as.character(round(auc_permuted, 3))))
  dev.off()

  list(permuted_confusionMatrix <- pred_permuted$cm)
}

rf.mtry_fit <- function(mtry, data, seed, trees){
  set.seed(seed)
  fit <- randomForest::randomForest(x = ( dplyr::select(data, -phenotype)),
                      y = data$phenotype,
                      data = data,
                      ntree = trees,
                      mtry = mtry,
                      importance = TRUE)
  fit[["err.rate"]][nrow(fit[["err.rate"]]),1]
}

rf.predict <- function (model, data){
  data$phenotype <- as.factor(data$phenotype)
  pred <- predict(model, newdata = data)
  cm <- confusionMatrix(data = pred, reference = data$phenotype)
  list(pred = pred, cm = cm)
}
