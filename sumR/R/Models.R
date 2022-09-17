#' @title Get names of models in Experiment
#' @param exp SummarizedExperiment object
#' @export
models <- function(exp){
  if (!validateExperiment(exp)) return(NULL)

  if (!"model" %in% names(metadata(exp))) return(NULL)
  return(names(metadata(exp)$model))
}

#' @title Get model in Experiment
#' @param exp SummarizedExperiment object
#' @param modelName Name of the model to retrieve
#' @export
model <- function(exp, modelName = 1){
  if (!validateExperiment(exp)) return(NULL)

  metadata(exp)$model[[modelName]]
}

#' @title Set value in model in Experiment
#' @param exp SummarizedExperiment object
#' @param modelName Name of the model to adjust
#' @param value Value to set in the model
#' @export
`model<-` <- function(exp, modelName, value){
  if (!validateExperiment(exp)) return(NULL)

  if (!"model" %in% names(metadata(exp))){
    metadata(exp)$model <- list()
  }
  metadata(exp)$model[[modelName]] <- value
  exp
}

#' @title Model using RandomForest model
#' @param exp SummarizedExperiment object
#' @param classifiers Column name in the colData that represents the classifiers
#' @param assay Which assay should be used to create a model
#' @param cv Number of folds to be used for cross validation
#' @param ratio Train-test split ratio
#' @param ... Any other parameters for preProcessing
#' @importFrom caret train createDataPartition trainControl varImp confusionMatrix
#' @importFrom stats predict
#' @export
generateModel <- function(exp, modelName, classifiers = metadata(exp)$phenotype,
                          assay = 1, folds = 5, ratio = 0.632, seed = NULL, ...){
  if (!validateExperiment(exp)) return(NULL)

  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  if (!is.null(seed)) set.seed(seed)

  data <- assay(exp, assay)

  trainIndex <- as.vector(createDataPartition(
    y = as.factor(exp[[classifiers]]), p = ratio, list = FALSE, times = 1
  ))

  modelList <- list()
  modelList$train <- colnames(exp)[trainIndex]
  modelList$test <- colnames(exp)[-trainIndex]

  control <- trainControl(method = "cv", number = folds,
                          savePredictions = "all",
                          preProcOptions = c("center", "scale"))

  modelList$model <- train(x = t(data[, modelList$train]),
                           y = as.factor(exp[, modelList$train][[classifiers]]),
                           method = modelName, trControl = control, ...)

  modelList$prediction <- predict(modelList$model,
                                  newdata = t(data[, modelList$test]),
                                  type = "prob")

  modelList$varImp <- varImp(modelList$model)

  x <- modelList$prediction
  pred <- factor(ifelse(x[,1] > 0.5, colnames(x)[1], colnames(x)[2]))
  obs <- factor(exp[, modelList$test][[metadata(exp)$phenotype]])
  modelList$confMatrix <- confusionMatrix(data = pred, reference = obs)

  model(exp, modelName) <- modelList
  exp
}

#' @title Variable Importance of models
#' @param exp SummarizedExperiment with model(s)
#' @param modelName name of the model
#' @export
varImportance <- function(exp, modelName = 1){
  if (!validateExperiment(exp)) return(NULL)

  df <- as.data.frame(model(exp, modelName)$varImp$importance)
  df$Compound <- rownames(df)
  df <- df[order(df$Overall, decreasing = TRUE), c("Compound", "Overall")]
  rownames(df) <- 1:nrow(df)
  df
}
