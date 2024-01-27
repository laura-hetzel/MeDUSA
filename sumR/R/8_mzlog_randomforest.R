# *** RandomForest Correlation -----------------------------------------------------
#' Find Correlation data within a mzLog_obj
#'
#' @returns Transposed mzlog_obj of correlation_data
#'
#' @export
mzlog_rf_correlation <- function(input_mzlog_obj, correlation_cutoff = 0.75){
  print("INFO: cor() is time and resource heavy on large data sets")
  data <- data.frame(t(dplyr::select(input_mzlog_obj,-mz)))
  high_cor_data <- caret::findCorrelation(cor(data), cutoff = correlation_cutoff)
  print(paste("INFO:mzlog_prep: ", length(high_cor_data)," of ", ncol(data) ," MZs are highly correlated ", sep=""))
  data <- data.frame(data[,-high_cor_data])
  colnames(data) <- input_mzlog_obj[-high_cor_data,"mz"]
  data <- cbind("sample_name"=rownames(data),data)
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
mzlog_rf_select <- function(correlation_data, metadata, attribute = "phenotype", feat_size_seq = seq(50,1000, by=50), plot = T) {
  data <- dplyr::left_join(correlation_data, metadata[c("sample_name", attribute)], by="sample_name")
  att <- as.factor(data[[attribute]])
  control <- caret::rfeControl( functions = caret::rfFuncs,
                                method = "repeatedcv",
                                repeats = 5,
                                number = 10)

  feat_select <- caret::rfe(dplyr::select(data, -"sample_name", -attribute),
                            att,
                            rfeControl = control,
                            sizes = feat_size_seq )
                            
  if(plot){
    ggplot(data = feat_select, metric = "Accuracy") + theme_bw()
    local.save_plot(paste("RandomForest Accuracy",local.mz_polarity_guesser(correlation_data),sep="-"))
    ggplot(data = feat_select, metric = "Kappa") + theme_bw()
    local.save_plot(paste("RandomForest Kappa",local.mz_polarity_guesser(correlation_data),sep="-"))
  }
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
#' @param attribute ##TODO Readd; currently only supports "phenotype" \cr
#'   String: which metadata attribute to compare 
#' @param seeds \cr
#'   List: which seeds to use. Also how many runs to do
#' @param cores \cr
#'   Int: can I haz multithreading
#'
#'
#' @export
mzlog_rf <- function(mzlog_obj,  metadata, rfe_obj = NULL, rf_trees = NULL, attribute = "phenotype", master_seeds = c(42,24),
                          seeds = c(201705,666,314159,1.05457, 998001), ratio = 0.8, cores = 4){
  if (! is.null(rfe_obj)){
    mzlog_obj <- mzlog_obj[ mzlog_obj$mz %in% as.numeric(caret::predictors(rfe_obj)), ]
  }
  if (is.null(rf_trees)){
    if( is.null(rfe_obj)){
      stop("ERROR: mzlog_rf_train: both rfe_obj & rf_trees are null")
    }
    rf_trees <- rfe_obj$bestSubset
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
  
  master_split <- .seed_splitter(data_t,attribute,master_seed[1],ratio)
  
  cl <- local.export_thread_env(cores, environment(mzlog_rf_train))
  tryCatch({
    out <- pbapply::pblapply(seeds, cl=cl, .run_train)
    out <- do.call(rbind,out)
    #plot ImpVars
  }, finally={
    local.kill_threads(cl)
  })
  out
  
  ##Validate
  train <- .transpose_massager(master_split$train, attribute)
  test <- .transpose_massager(master_split$test, attribute)


  set.seed(master_seed[2])
  model <- randomForest(x = train_neg[-1],
                        y = train_neg[[attribute]],
                        data = train_neg,
                        ntree = rf_tree,
                        mtry = mtry,
                        importance = TRUE, 
                            proximity = TRUE)
}

.transpose_masssager <- function(data, attribute){
  rand <- sample(1:nrow(data),nrow(data))
  data[[attribute]] <- paste(data[[attribute]], rand, sep="_")
  rownames(data) <- NULL
  data <- column_to_rownames(data, attribute)
  
  data <- as.data.frame(t(data))
  data <- rownames_to_column(data, "mz")
  data$mz <- as.numeric(gsub("X","",data$mz))
  data
}

.seed_splitter <- function(data,attribute,seed, ratio){
  set.seed(seed)
  split <- caTools::sample.split(data[[attribute]], SplitRatio = ratio)
  tr <- subset(data, split == T)
  te <- subset(data, split == F)
  list(train = tr, test = te)
}

.run_train <- function( seed, master_split, ratio, rf_trees, attribute, mtry, tunegrid, control){
  #Split
  split <- .seed_splitter(master_split$train, attribute, seed, ratio)
  #TODO Dumb workaround to make "as.factor()" work: (it doesn't like parameters)
  tmp <- split$train
  colnames(tmp)[colnames(tmp) == attribute] <- "llama"  
  #Train
  rf_fit <- caret::train( as.factor(llama) ~.,
                   data = tmp,
                   method = 'rf',
                   tuneGrid = tunegrid,
                   trControl = control,
                   ntree = rf_trees ,
                   na.action = na.exclude)
  
  rf_pred <- predict(rf_fit, split$test)
  acc <- caret::confusionMatrix(rf_pred, as.factor(split$test[[attribute]]))$overall
  #TODO
  # polarity guesser
  # AOC
  # ROC
  # Accuracty
  imp <- caret::varImp(rf_fit)$importance
  imp <- rownames_to_column(imp, "mz")
  imp$mz <- as.numeric(gsub("X", "", imp$mz))
  imp <- head(imp[order(-imp$Overall),], n=100)
  rownames(imp) <- NULL
  out <- imp
}
