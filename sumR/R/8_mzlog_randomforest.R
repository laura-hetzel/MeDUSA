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
  metadata <- local.meta_polarity_fixer(correlation_data, metadata)
  data <- dplyr::left_join(correlation_data, meta[c("sample_name", attribute)], by="sample_name")
  att <- as.factor(data[[attribute]])
  control <- caret::rfeControl( functions = caret::rfFuncs,
                                method = "repeatedcv",
                                repeats = 5,
                                number = 10)

  feat_select <- caret::rfe(dplyr::select(data, -sample_name, -attribute),
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


# *** RandomForest model -----------------------------------------------------
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
#' @returns list(model, test, train)
#'
#' @export
mzlog_rf <- function(mzlog_obj,  metadata, rfe_obj = NULL, rf_trees = NULL, attribute = "phenotype", master_seeds = c(42,1701),
                          seeds = c(201705, 623.1912, 314159, 1.05457, 10191), ratio = 0.8, cores = 4, plot = T){
  if (! is.null(rfe_obj)){
    mzlog_obj <- mzlog_obj[ mzlog_obj$mz %in% as.numeric(caret::predictors(rfe_obj)), ]
  }
  if (is.null(rf_trees)){
    if( is.null(rfe_obj)){
      stop("ERROR: mzlog_rf_train: both rfe_obj & rf_trees are null")
    }
    rf_trees <- rfe_obj$bestSubset
  }
  metadata <- local.meta_polarity_fixer(mzlog_obj, metadata)
  pol <- local.mz_polarity_guesser(mzlog_obj)
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

  master_split <- rf.seed_splitter(data_t,attribute,master_seeds[1],ratio)

  #TODO, better manage variables among threads
  rf_vars <- list("master_split" = master_split, "ratio" = ratio,
                  "rf_trees" = rf_trees, "attribute" = attribute, "pol" = pol,
                  "mtry"= mtry, "tunegrid" = tunegrid, "control" = control)

  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    out <- pbapply::pblapply(seeds, cl=cl, rf.run_train, rf_vars = rf_vars)
    out <- do.call(rbind,out)
    #plot ImpVars
  }, finally={
    local.kill_threads(cl)
  })
  out

  ##Validate
  train <- rf.transpose_massager(master_split$train, attribute)
  test <- rf.transpose_massager(master_split$test, attribute)

  set.seed(master_seeds[2])
  model <- randomForest(x = master_split$train[-1],
                        y = master_split$train[[attribute]],
                        data = master_split$train,
                        ntree = rf_trees,
                        mtry = mtry,
                        importance = TRUE,
                        proximity = TRUE)

  distance_matrix <- as.dist(1 - model$proximity)
  mds <- cmdscale(distance_matrix, eig = TRUE, x.ret = TRUE)
  mds_var <- round(mds$eig/sum(mds$eig)*100, 1)
  mds_values <- mds$points
  mds_data <- data.frame(Sample = rownames(mds_values),
                             x = mds_values[,1],
                             y = mds_values[,2],
                             Phenotype = master_split$train$phenotype)
  if(plot){
    plot(model)
    local.save_plot(paste0("rf_model_",pol))

    varImpPlot(model, col = "Blue", pch = 2, main = paste("Important Variables ",pol))
    local.save_plot(paste0("rf_varImp_",pol))
    pheno <- as.vector(unique(master_split$train$phenotype))
    ggplot(mds_xxdata, aes(x = x, y = y, label = Sample)) +
      geom_point(aes(color = Phenotype), size = 3) +
      theme_classic() +
      theme(legend.position = "bottom",) +
      theme(text = element_text(size = 18)) +
      xlab(paste("MDS1 - ", mds_var_neg[1], "%", sep = "")) +
      ylab(paste("MDS2 - ", mds_var_neg[2], "%", sep = "")) +
      scale_color_manual(breaks = pheno,
                         values = c("red4", "deepskyblue2")) +
      ggtitle(paste0("MDS plot ", pol," using random forest proximities"))
      local.save_plot(paste0("rf_mds_",pol))

  }

  list( model = model, test = master_split$test, train = master_split$train )
}


rf.transpose_massager <- function(data, attribute){
  rand <- sample(1:nrow(data),nrow(data))
  data[[attribute]] <- paste(data[[attribute]], rand, sep="_")
  rownames(data) <- NULL
  data <- tibble::column_to_rownames(data, attribute)

  data <- as.data.frame(t(data))
  data <- tibble::rownames_to_column(data, "mz")
  data$mz <- as.numeric(gsub("X","",data$mz))
  data
}

rf.seed_splitter <- function(data, attribute, seed, ratio){
  set.seed(seed)
  split <- caTools::sample.split(data[[attribute]], SplitRatio = ratio)
  tr <- subset(data, split == T)
  te <- subset(data, split == F)
  list(train = tr, test = te)
}

rf.run_train <- function( seed, rf_vars){
  #Split
  split <- rf.seed_splitter(rf_vars$master_split$train, rf_vars[[attribute]], seed, rf_vars$ratio)
  #TODO Dumb workaround to make "as.factor()" work: (it doesn't like parameters)
  tmp <- split$train
  colnames(tmp)[colnames(tmp) == rf_vars$attribute] <- "factor"
  #Train
  rf_fit <- caret::train( as.factor(factor) ~.,
                   data = tmp,
                   method = 'rf',
                   tuneGrid = rf_vars$tunegrid,
                   trControl = rf_vars$control,
                   ntree = rf_vars$rf_trees ,
                   na.action = na.exclude)

  rf_pred <- predict(rf_fit, split$train$phenotype )
  acc <- caret::confusionMatrix(rf_pred, as.factor(split$test[[rf_vars$attribute]]))$overall
  roc <- pROC::roc(tmp$phenotype, rf_pred)
  auc <- pROC::auc(roc)
  plot(roc, col = "Red", main = paste0("Fold_",pol,"_", seed),
       sub = paste0("Acc:",acc," AUC:", as.character(round(auc, 3))))
  local.save_plot(paste0("rf_train_",seed))

  imp <- caret::varImp(rf_fit)$importance
  imp <- tibble::rownames_to_column(imp, "mz")
  imp$mz <- as.numeric(gsub("X", "", imp$mz))
  imp <- head(imp[order(-imp$Overall),], n=100)
  rownames(imp) <- NULL
  out <- imp
}
