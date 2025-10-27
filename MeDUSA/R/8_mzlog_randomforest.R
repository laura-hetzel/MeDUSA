# *** RandomForest Correlation -----------------------------------------------------
#' MzLog: RandomrForest generate correlation data
#'
#' @description
#' Random Forest is a robust tool for identifying the features that contribute
#' to correct phenotype prediction. For optimal performance of the model, it is
#' recommended to remove the highly correlated features before training and
#' testing the model. The mzlog_rf_correlation function utilizes the
#' findCorrelation function of the Caret package to isolate and remove the
#' highly correlated features. The new data set with these feature removed will
#' be used for the remaining functions of random forest.
#'
#' @param input_mzlog_obj \cr
#'   DataFrame : MzLog object
#' @param correlation_cutoff \cr
#'   Float: Decimal percentage. Higher cutoff = less correlation
#' @examples
#'   mzlog_rf_correlation(input_mzlog_obj) : run standard correlation
#'   mzlog_rf_correlation(input_mzlog_obj, correlation_cutoff = 0.8) : run correlation with custom cutoff percentage
#' @return correlation_data: a Transposed mzlog_obj of non-correlated_data
#' @export
mzlog_rf_correlation <- function(input_mzlog_obj, correlation_cutoff = 0.75){
  data <- data.frame(t(dplyr::select(input_mzlog_obj,-mz)))
  print("INFO: cor() can be time and resource heavy on large data sets; 20k mzs can take 18gb memory")
  cor <- cor(data)
  high_cor_data <- caret::findCorrelation(cor , cutoff = correlation_cutoff)
  if ( ncol(data) - length(high_cor_data) < 2 ) {
    stop("ERROR:MeDUSA::mzlog_rf_correlation: Would return less than 2 non-correlated samples; try increasing 'correlation_cutoff'")
  }
  print(paste("INFO:mzlog_prep: ", length(high_cor_data)," of ", ncol(data) ," MZs are highly correlated ", sep=""))
  data <- data.frame(data[,-high_cor_data])
  colnames(data) <- input_mzlog_obj[-high_cor_data,"mz"]
  data <- cbind("sample_name"=rownames(data),data)
}

# *** RandomForest Select -----------------------------------------------------
#' correlation_data: RandomrForest rfe_select
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
#' processing.
#'
#' @param correlation_data \cr
#'   DataFrame : from rf_correlation (is a transposed mz_obj)
#' @param metadata \cr
#'   DataFrame: metadata object
#' @param feat_size_seq \cr
#'   Sequence to find optimal "number_of_variables"
#' @param plot \cr
#'   Boolean   : To plot or not to plot.
#' @examples
#'   mzlog_rf_select(correlation_data, metadata) : run standard rf_select
#'   mzlog_rf_select(correlation_data, metadata, seq(100,1000), by=100)) : run rf_seelct with custom feature sizes
#' @return caret::rfe object
#' @export
mzlog_rf_select <- function(correlation_data, metadata, feat_size_seq = seq(50,1000, by=50), plot = T) {
  metadata <- local.meta_polarity_fixer(correlation_data, metadata)
  data <- dplyr::left_join(correlation_data, metadata[c("sample_name", "phenotype")], by="sample_name")
  att <- as.factor(data$phenotype)
  control <- caret::rfeControl( functions = caret::rfFuncs,
                                method = "repeatedcv",
                                repeats = 5,
                                number = 10)

  print("INFO:MeDUSA::mzlog_rf_select: caret::rfe() is time and resource heavy on large data sets")
  feat_select <- caret::rfe(dplyr::select(data, -sample_name, -phenotype),
                            att,
                            rfeControl = control,
                            sizes = feat_size_seq )

  if(plot){
    ggplot(data = feat_select, metric = "Accuracy") + theme_bw()
    local.save_plot(paste0("RandomForest Accuracy-",local.mz_polarity_guesser(correlation_data)))
    ggplot(data = feat_select, metric = "Kappa") + theme_bw()
    local.save_plot(paste0("RandomForest Kappa-",local.mz_polarity_guesser(correlation_data)))
  }
  feat_select
}

# *** RandomForest model -----------------------------------------------------
#' MzLog & rfe_select: Run random forest
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
#' fundamental the feature was in predicting the phenotype. The accuracy of the
#' model is also output.
#'
#' @param mzlog_obj \cr
#'   DataFrame : mzlog : expects "phenotype"
#' @param metadata \cr
#'   DataFrame: metadata object
#' @param rfe_obj \cr
#'   carat::rfe_obj (see rf_select)
#' @param feat_size_seq \cr
#'   Sequence to find optimal "number_of_variables"
#' @param master_seeds \cr
#'   List: Two seeds to use for training
#' @param seeds \cr
#'   List: How many training runs to do, and which seeds to use.
#' @param ratio \cr
#'   List: Ratio of training vs test data
#' @param cores \cr
#'   Int: Can I haz multi-threading
#' @param plot \cr
#'   Bool: To plot or not to plot?
#' @examples
#'   mzlog_rf(input_mzlog_obj, metadata, rfe_obj = rf_select_ouput) : run recommeneded random forest
#'   mzlog_rf(input_mzlog_obj, metadata, seeds = c(1.23, 45.6, 789)) : run for only 3 runs
#' @return list(model, test, train, imp_mz)
#' @export
mzlog_rf <- function(mzlog_obj,  metadata, rfe_obj = NULL, rf_trees = NULL, master_seeds = c(42,1701),
                          seeds = c(201705, 623.1912, 314159, 1.05457, 10191), ratio = 0.8, cores = 4, plot = T){
  if (! is.null(rfe_obj)){
    mzlog_obj <- mzlog_obj[ mzlog_obj$mz %in% as.numeric(caret::predictors(rfe_obj)), ]
  }
  if (is.null(rf_trees)){
    if( is.null(rfe_obj)){
      stop("ERROR:MeDUSA::mzlog_rf: both rfe_obj & rf_trees are null")
    }
    rf_trees <- rfe_obj$bestSubset
  }
  metadata <- local.meta_polarity_fixer(mzlog_obj, metadata)
  pol <- local.mz_polarity_guesser(mzlog_obj)
  rownames(mzlog_obj) <- mzlog_obj$mz
  data_t <- data.frame(t(dplyr::select(mzlog_obj,-mz)))
  data_t <- tibble::rownames_to_column(data_t,"sample_name")
  data_t <- dplyr::left_join(data_t, metadata[c("sample_name", "phenotype")])
  data_t <- dplyr::select(data_t,-"sample_name")
  data_t$phenotype <- as.factor(data_t$phenotype)

  #Train settings
  mtry <- c(sqrt(ncol(data_t)))
  tunegrid <- expand.grid(.mtry=mtry)
  control <- trainControl(method ='repeatedcv',
                          number = 10,
                          repeats = 4,
                          search = 'grid',
                          allowParallel = TRUE)
  master_split <- rf.seed_splitter(data_t,master_seeds[1],ratio)

  #TODO, better manage variables among threads
  rf_vars <- list("master_split" = master_split,
                  "rf_trees" = rf_trees, "pol" = pol, "ratio" = ratio,
                  "mtry"= mtry, "tunegrid" = tunegrid, "control" = control)

  tryCatch({
    cl <- local.export_thread_env(cores, environment())
    out <- pbapply::pblapply(seeds, cl=cl, rf.run_train, rf_vars = rf_vars)
    out <- do.call(rbind,out)
  }, finally={
    local.kill_threads(cl)
  })
  out <- out[order(-out$Overall),]
  out = out[!duplicated(out$mz),]
  ##Validate
  train <- rf.transpose_massager(master_split$train)
  test <- rf.transpose_massager(master_split$test)

  set.seed(master_seeds[2])
  model <- randomForest(x = select(master_split$train, -phenotype),
                        y = master_split$train$phenotype,
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
    tmp <- as.data.frame(model$err.rate)
    tmp$tree <- c(1:nrow(tmp))
    tmp <- tidyr::pivot_longer(tmp, !tree )
    ggpubr::ggline(tmp,x="tree",y="value", shape = "", colro = "name", line = "name", title=paste0("RandomForest Error: ",pol))
    local.save_plot(paste0("RF_error", pol))

    png(paste0(local.output_dir(),local.dir_sep(),"Important_Variables_",pol,".png"))
    varImpPlot(model, col = "Blue", pch = 2, main = paste0("Important Variables ", pol))
    dev.off()

    pheno <- as.vector(unique(master_split$train$phenotype))
    ggplot(mds_data, aes(x = x, y = y, label = Sample)) +
      geom_point(aes(color = Phenotype), size = 3) +
      theme_classic() +
      theme(legend.position = "bottom",) +
      theme(text = element_text(size = 18)) +
      xlab(paste0("MDS1 - ", mds_var[1], "%")) +
      ylab(paste0("MDS2 - ", mds_var[2], "%")) +
      scale_color_manual(breaks = pheno,
                         values = c("red4", "deepskyblue2")) +
      ggtitle(paste0("MDS plot ", pol," using random forest proximities"))
      local.save_plot(paste0("randomForest_mds_",pol))
  }

  list( model = model, test = master_split$test, train = master_split$train, imp_mz = out )
}

# *** RandomForest  -----------------------------------------------------------
#' MzLog: Run random forest
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
#' fundamental the feature was in predicting the phenotype. The accuracy of the
#' model is also output. The rf_validate function is used to compare the
#' predictions of the model to the actual phenotype of the sample. The accuracy
#' of the model is also output. To ensure the accuracy of the model is truly
#' based on biological phenotype, a permutation of the data is used. The
#' rf_permuted function randomly assigns a phenotype to each of the samples and
#' feeds this permuted data set through the model. If the accuracy is
#' significantly greater than 50%, the user should reconsider the model.
#' The mzlog_rf_magic function incorporates all of the mentioned functions of
#' random forest into one seamless function for the user.
#'
#' @param mzlog_obj \cr
#'   DataFrame : mzlog : expects "phenotype"
#' @param metadata \cr
#'   DataFrame: metadata object
#' @param polarity \cr
#'   Character: "pos" | "neg"
#' @examples
#'   mzlog_rf_magic(input_mzlog_obj, metadata, "pos") : run random forest for postive
#' @return List of important mz
#' @export
mzlog_rf_magic <- function(input_mzlog_obj, metadata, polarity){
  tryCatch({
    rf_cor <- mzlog_rf_correlation(input_mzlog_obj)
    gc()
    print("INFO:MeDUSA::mzlog_rf_magic: Correlation success")
  }, error = function(e) {
    print(e)
    stop("ERROR: in mzlog_rf_correlation")
  })
  tryCatch({
    rf_sel <- mzlog_rf_select(rf_cor, metadata)
    gc()
    print("INFO:MeDUSA::mzlog_rf_magic: Feature Select success")
  }, error = function(e) {
    print(e)
    stop("ERROR:MeDUSA::mzlog_rf_magic: in mzlog_rf_select")
  })
  tryCatch({
    rf_obj <- mzlog_rf(input_mzlog_obj, metadata, rfe_obj = rf_sel )
    gc()
    print("INFO:MeDUSA::mzlog_rf_magic: RandomForest success")
  }, error = function(e) {
    print(e)
    stop("ERROR:MeDUSA::mzlog_rf_magic: in mzlog_rf")
  })
  tryCatch({
    rf_validate <- rf_validate(rf_obj, trees = rf_sel$bestSubset * 1.1)
    gc()
    print("INFO:MeDUSA::mzlog_rf_magic: RF validation success")
  }, error = function(e) {
    print(e)
    stop("ERROR:MeDUSA::mzlog_rf_magic: in rf_validate")
  })
  tryCatch({
    rf_permuted <- rf_permuted(rf_obj, polarity)
    gc()
    print("INFO:MeDUSA::mzlog_rf_magic: RF permuted success")
  }, error = function(e) {
    print(e)
    stop("ERROR:MeDUSA::mzlog_rf_magic: in rf_permuted")
  })
  rf_obj$imp_mz
}

rf.transpose_massager <- function(data){
  rand <- sample(1:nrow(data),nrow(data))
  data$phenotype <- paste(data$phenotype, rand, sep="_")
  rownames(data) <- NULL
  data <- tibble::column_to_rownames(data, "phenotype")

  data <- as.data.frame(t(data))
  data <- tibble::rownames_to_column(data, "mz")
  data$mz <- as.numeric(gsub("X","",data$mz))
  data
}

rf.seed_splitter <- function(data, seed, ratio){
  set.seed(seed)
  split <- caTools::sample.split(data$phenotype, SplitRatio = ratio)
  tr <- subset(data, split == T)
  te <- subset(data, split == F)
  list(train = tr, test = te)
}

rf.run_train <- function( seed, rf_vars){
  #Split
  split <- rf.seed_splitter(rf_vars$master_split$train, seed, rf_vars$ratio)
  #Train
  rf_fit <- caret::train( as.factor(phenotype) ~.,
                   data =  split$train,
                   method = 'rf',
                   tuneGrid = rf_vars$tunegrid,
                   trControl = rf_vars$control,
                   ntree = rf_vars$rf_trees ,
                   na.action = na.exclude)
  rf_pred <- predict(rf_fit, split$test )
  acc <- caret::confusionMatrix(rf_pred, as.factor(split$test$phenotype))
  roc <- pROC::roc(as.numeric(split$test$phenotype), as.numeric(rf_pred))
  auc <- pROC::auc(roc)

  rf_train_data <- list('acc'=acc, 'roc'=roc, 'auc'=auc)
  save(rf_train_data, file = paste0(local.output_dir(), local.dir_sep(), "RF_train", pol, "_", seed,".Rdata"))

  pROC::ggroc(roc, color='red') +
      ggtitle( label = paste0("ROC: ", pol), subtitle = paste0("Seed:", seeds[1]))
  local.save_plot(paste("RF_ROC",pol, seed, sep="_"))

  imp <- caret::varImp(rf_fit)$importance
  imp <- tibble::rownames_to_column(imp, "mz")
  imp$mz <- as.numeric(gsub("X", "", imp$mz))
  imp <- head(imp[order(-imp$Overall),], n=100)
  rownames(imp) <- NULL
  out <- imp
}
