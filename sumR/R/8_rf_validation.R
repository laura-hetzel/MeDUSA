
# *** RandomForest Validate-----------------------------------------------------
#' Train you model data within a mzLog_obj
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
rf_validate <- function(rf_obj, mtry_range = c(1:200), mtry_seed = 1984, cores = 2){
  pred_train <- rf.predict(rf_obj$model, rf_obj$train)
  if( sum(!(pred_train$pred == rf_obj$train$phenotype)) ){
    warning("sumR:: rf_validate: Train prediction did not match train model")
  }
  pred_test <- rf.predict(rf_obj$model, rf_obj$test)
  results <- as.data.frame(cbind("actual" = rf_obj$test$phenotype,
                                 "prediction" = pred_test$pred))

  ml <- data.frame(mtry = mtry_range)
  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    ml$err_rate <- pbapply::pbapply(ml, 1, rf.mtry_fit, cl = cl,
                                   data = rf_obj$train, seed = mtry_seed)
    best_mtry <- ml[ml$err_rate == min(ml$err_rate),]
    print(paste0("INFO: sumR::rf_validate: Lowest error:", best_mtry$err_rate,
                " At mtry:", best_mtry$mtry ))
  }, finally={
   local.kill_threads(cl)
  })

  list(train_confusionMatrix <- pred_train$cm, test_confusionMatrix <- pred_test$cm, mtry_error <- ml)
}


# *** RandomForest Permuted -----------------------------------------------------
#' Train you model data within a mzLog_obj
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
rf_permuted <- function(rf_obj, pol, seed = 105.9, plot = T){
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

  plot(roc_permuted, col = "Red", main = paste0(pol,"_RF_ROC Permuted"),
       sub = paste0("Acc:",pred_permuted$cm$overall["Accuracy"]," AUC:", as.character(round(auc_permuted, 3))))
  local.save_plot(paste0("randomForest_ROC_Permuted_",pol))


  list(permuted_confusionMatrix <- pred_permuted$cm)
}

rf.mtry_fit <- function(mtry, data, seed){
  set.seed(seed)
  fit <- randomForest::randomForest(x = ( dplyr::select(data, -phenotype)),
                      y = data$phenotype,
                      data = data,
                      ntree = 500,
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
