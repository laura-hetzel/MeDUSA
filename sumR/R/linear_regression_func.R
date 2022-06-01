#' @title Model glmnet
#' @export
modelGLM <- function(exp, classifiers = metadata(exp)$phenotype, assay = 1,
                        type.measure = "deviance", lambda = c("lambda.min", "lambda.1se"),
                        cv = 5, ratio = 0.8, seed = NULL, alphas = seq(0, 1, 0.1)){
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  if (!is.null(seed)) set.seed(seed)

  data <- assay(exp, assay)
  samples <- as.factor(exp[[classifiers]])
  split <- caTools::sample.split(samples, SplitRatio = ratio)
  model(exp, "glmnet")$train <- colnames(exp)[which(split)]
  model(exp, "glmnet")$test <- colnames(exp)[-which(split)]


  training <- t(data[, model(exp, "glmnet")$train])
  test <- t(data[, model(exp, "glmnet")$test])

  folds <- caret::createFolds(y = model(exp, "glmnet")$train,
                              k = cv, list = F)
  alpha <- alphas[which.min(vapply(alphas, function(a){
    tryCatch({
      return(
        cv.glmnet(
          y = exp[, model(exp, "glmnet")$train][[classifiers]],
          x = as.matrix(scale(training)), foldid = folds,
          family = "binomial", alpha = a, type.measure = type.measure,
          nfolds = cv
        )[[lambda[1]]]
      )
    }, error = function(x) return(Inf))

  }, double(1)))]
  m <- cv.glmnet(
    y = exp[, model(exp, "glmnet")$train][[classifiers]],
    x = as.matrix(scale(training)), foldid = folds,
    family = "binomial", alpha = alpha, type.measure = type.measure,
    nfolds = cv
  )
  model(exp, "glmnet")$model <- m

  model(exp, "glmnet")$predict <- predict(m, s = m[[lambda[1]]],
                                          newx = as.matrix(scale(test)), type = "class")
  exp
  # do.call(rbind, lapply(0:10, function(i){
  #   m <- cv.glmnet(
  #     y = exp[, model(exp, "glmnet")$train][[classifiers]],
  #     x = as.matrix(scale(training)), foldid = folds,
  #     family = "binomial", alpha = i / 10, type.measure = type.measure,
  #     nfolds = cv
  #   )
  #
  #
  #   predict <- predict(m, s = m[[lambda[1]]],
  #                      newx = as.matrix(scale(test)),
  #                      type = type.measure
  #   )
  #
  #
  #   cf <- confusionMatrix(table(data = predict, reference = exp[, model(exp, "glmnet")$test][[classifiers]]))
  #   data.frame(alpha = i / 10, cf_acc = cf$overall$Accuracy)
  # }))
}

## choosing alpha value for elastic net
#' @title Alpha value choice
#' @description This function choose the alpha value for elastic net regression with the desired accuracy
#' @param training training set
#' @param test test set
#' @param seed global seed for reproducible results
#' @param type.measure loss to use for cross-validation. type.measure="class" ,"auc"
#' @param lambda choose between "lambda.min" or "lambda.1se"
#' @importFrom glmnet cv.glmnet
#' @importFrom caret confusionMatrix
#' @importFrom dplyr select
glmnet_cv_alphachoice <- function(training, test, type.measure, seed = NULL, lambda = c("lambda.min", "lambda.1se")) {
  #list_fits <- list()
  if (!is.null(seed)) set.seed(seed)

  do.call(rbind, lapply(0:10, function(i){
    model <- cv.glmnet(
      y = training$samples,
      x = as.matrix(scale(training)),
      family = "binomial", alpha = i / 10, type.measure = type.measure
    )
    predict <- predict(model, s = model[[lambda[1]]],
                       newx = as.matrix(scale(test)),
                       type = type.measure
    )
    cf <- confusionMatrix(table(data = predict, reference = test$samples))
    data.frame(alpha = i / 10, cf_acc = cf$overall$Accuracy)
  }))
}


#
#' @title final glmnet model
#' @description This function makes the final model of glmnet using the desired alpha value
#' @param training training set
#' @param test test set
#' @param seed global seed for reproducible results
#' @param type.measure loss to use for cross-validation. type.measure="class" ,"auc"
#' @param nfolds number of folds for cross validation
#' @importFrom glmnet cv.glmnet
#' @importFrom dplyr select
final_glmnet <- function(training, test, alpha, type.measure, seed = NULL, nfolds = 10) {
  if (!is.null(seed)) set.seed(seed)
  cv.glmnet(
    y = training$samples,
    x = as.matrix(scale(training %>% dplyr::select(-samples))),
    family = "binomial", alpha = alpha, type.measure = type.measure,
    keep = T, nfolds = nfolds
  )
}


#' @title final glmnet model plot
#' @param final_fit the cv.glmnet model
final_glmnet_plot <- function(fit_final) {
  plot(fit_final)
}


#' @title final glmnet model assess
#' @param final_fit the cv.glmnet model
#' @param test test set
#' @param lambda the choice of lambda . default is "lambda.min"
#' @importFrom dplyr select
#' @importFrom glmnet assess.glmnet
final_glmnet_assess <- function(fit_final, test, lambda = "lambda.min") {
  assess.glmnet(fit_final, newx = as.matrix(scale((test %>% dplyr::select(-samples)))),
                newy = test$samples, s = lambda)
}



#' @title final glmnet model confusion matrix
#' @param final_fit the cv.glmnet model
#' @param test test set
#' @param lambda the choice of lambda . default is "lambda.min"
#' @importFrom caret confusionMatrix
#' @importFrom dplyr select
final_glmnet_confusionmatrix <- function(fit_final, test, lambda = "lambda.min") {
  predict <- predict(fit_final, s = lambda,
                     newx = as.matrix(scale((test %>% dplyr::select(-samples)))),
                     type = "class")
  confusionMatrix(table(data = predict, reference = test$samples))
}



#' @title final glmnet model confusion matrix
#' @param final_fit the cv.glmnet model
#' @param master_df the master df of the raw or scaled intensitites that have m/z column
#' @param lambda the choice of lambda . default is "lambda.min"
#' @importFrom dplyr filter_all
#' @importFrom dplyr any_vars
final_glmnet_features <- function(fit_final, master_df, lambda = "lambda.min") {
  k <- coef(fit_final, s = lambda)
  k2 <- k[which(k != 0), ]
  k2
  k3 <- as.data.frame(names(k2))
  k3 <- k3[-1, ]
  imp_var <- master_df %>% filter_all(any_vars(. %in% k3)) ## obtaining the intensities of the all imp mz from the master df
  return(list("selected_feature" = imp_var, "coefficients" = k))
}




#
#' @title roc curve for auc
#' @description This function makes ROC curve for AUC measure of cv.glmnet
#' @param training training set
#' @param final_fit the cv.glmnet model with "auc" as type.measure
#' @importFrom glmnet roc.glmnet
roc_binomial <- function(final_fit, training) {
  rocs <- roc.glmnet(final_fit$fit.preval, newy = training$samples)
  best <- final_fit$index["min", ]
  plot(rocs[[best]], type = "l")
  invisible(sapply(rocs, lines, col = "grey"))
  lines(rocs[[best]], lwd = 2, col = "red")
}
