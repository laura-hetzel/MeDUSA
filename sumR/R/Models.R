#' @title Get names of models in Experiment
#' @description After modelling, the results are stored in the [metadata] slot
#' as a list with name `model`. However, multiple models can be made and stored
#' this way. This function returns the names of models that have been created.
#' @details The `model` name is a reserved name in sumR, as all models are
#' stored in a list with the name `model`. Created models and their results like
#' cross-validation and variable importances are stored in a sub-list of
#' `model`. These can be accessed with the [model] function, which requires
#' a model name or index. This function aids in this process by returning the
#' names of all models that have been created with [generateModel].
#' @returns Character vector of the names of models created with
#' [generateModel].
#' @param exp SummarizedExperiment object
#' @export
#' @examples
models <- function(exp){
  if (!validateExperiment(exp)) return(NULL)

  if (!"model" %in% names(metadata(exp))) return(NULL)
  return(names(metadata(exp)$model))
}

#' @title Get model in Experiment
#' @description After modelling with [generateModel], the results are stored in
#' the [metadata] slot as a list with name `model`. However, multiple models
#' can be made and stored this way. This function returns the model information
#' of the given model, either by index or name of the model used.
#' @inherit generateModel details
#' @returns List of various slots to inspect the model. Slots include:
#' * __train__: Samples that are used as a train set
#' * __test__: Samples that are used as a test set. These won't be taken into
#' account during cross-validation to prevent bias. The results from testing the
#' model on these samples is stored in the `prediction` and `confMatrix` slot.
#' * __control__: Settings that were used to train the model, including for
#' cross-validation.
#' * __model__: The best model that was picked during cross-validation.
#' * __prediction__: Result of the [predict] function on the test dataset.
#' * __varImp__: The variable importance of the model during training. Can be
#' used to assess what combination of variables are the best predictors for
#' the phenotype. These results can also be retrieved using the [varImportance]
#' function.
#' * __confMatrix__: Confusion matrix of the predictions of the model on the
#' testset. Can be used to assess the performance of the model.
#' @seealso [generateModel], [varImportance], [models]
#' @param exp SummarizedExperiment object obtained after post-processing.
#' @param modelName Name or index of the model to retrieve. Defaults to `1`
#' @export
#' @examples
model <- function(exp, modelName = 1){
  if (!validateExperiment(exp)) return(NULL)

  metadata(exp)$model[[modelName]]
}

#' @title Set or replace a model in the given experiment
#' @description
#' @details
#' @returns
#' @param exp SummarizedExperiment object
#' @param modelName Name of the model to adjust
#' @param value Value to set in the model
`model<-` <- function(exp, modelName, value){
  if (!validateExperiment(exp)) return(NULL)

  if (!"model" %in% names(metadata(exp))) {
    metadata(exp)$model <- list()
  }
  metadata(exp)$model[[modelName]] <- value
  exp
}

#' @title Model phenotypes using the caret library
#' @description This function trains and tests a given model from the caret
#' library to predict phenotypes. This can be used as a multivariate approach
#' to determine important differentiating peaks between phenotypes.
#' @details Modelling is a multi-variate approach of determining important
#' peaks to predict the phenotype(s). Models are created using the [caret]
#' package. First, a train-test split is made using the given ratio. Next,
#' train control is done using a cross-validation approach. Here, the data
#' is scaled and centered before training starts. The model used for training
#' can be picked from the caret library of classification models. The results
#' are stored in the model field of the SummarizedExperiment metadata slot. It
#' can be accessed by the [model] function.
#' @returns SummarizedExperiment with updated model field in the metadata slot.
#' @param exp SummarizedExperiment object
#' @param classifiers Column name in the colData that represents the classifiers
#' @param assay Which assay should be used to create a model
#' @param cv Number of folds to be used for cross validation
#' @param ratio Train-test split ratio
#' @param ... Any other parameters for preProcessing
#' @importFrom caret train createDataPartition trainControl varImp confusionMatrix
#' @importFrom stats predict
#' @seealso [Caret manual](https://topepo.github.io/caret/)
#' @export
#' @examples
generateModel <- function(exp, modelName, classifiers = metadata(exp)$phenotype,
                          assay = 1, folds = 5, ratio = 0.632, seed = NULL, ...){
  if (!validateExperiment(exp)) return(NULL)

  if (is.null(classifiers)) stop("Cannot perform test without phenotype")
  if (!is.null(seed)) set.seed(seed)

  data <- assay(exp, assay)

  trainIndex <- as.vector(createDataPartition(
    y = as.factor(exp[[classifiers]]), p = ratio, list = FALSE, times = 1
  ))

  modelList <- list()
  modelList$train <- colnames(exp)[trainIndex]
  modelList$test <- colnames(exp)[-trainIndex]

  modelList$control <- trainControl(method = "cv", number = folds,
                          savePredictions = "all",
                          preProcOptions = c("center", "scale"))

  modelList$model <- train(x = t(data[, modelList$train]),
                           y = as.factor(exp[, modelList$train][[classifiers]]),
                           method = modelName, trControl = modelList$control,
                           ...)

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

#' @title Assess the variable importance of models
#' @description This function retrieves that variable importances of peaks
#' that were used in a model. For each peak, the probability is shown in the
#' `Overall` column.
#' @details Variable importance is a metric that determines how important the
#' variable (compound) was for determining the difference between the given
#' phenotypes. A higher probability indicates a larger difference and a
#' stronger association with a phenotype.
#' @returns A dataframe with two columns: `Compound` and `Overall`, which
#' indicates the overall probability of a compound for a model
#' @param exp SummarizedExperiment with model(s) generated by [generateModel].
#' @param modelName Name or index of the model. Defaults to `1`, meaning the
#' first model that was generated.
#' @export
#' @examples
varImportance <- function(exp, modelName = 1){
  if (!validateExperiment(exp)) return(NULL)

  df <- as.data.frame(model(exp, modelName)$varImp$importance)
  df$Compound <- rownames(df)
  df <- df[order(df$Overall, decreasing = TRUE), c("Compound", "Overall")]
  rownames(df) <- 1:nrow(df)
  df
}
