## shapiro's test for normality check
## dataframe is transposed with mz values as columns

#' @title shapiro test 
#' @description shapiro test for normality check
#' @param dataframe transposed dataframe with m/z as columns 
#' @param threshold numerical value for wanted threshold . the default is 0.05
shapiro_results <- function(dataframe,threshold=0.05){
  Test_results <- apply(dataframe,2, function(x) shapiro.test(as.numeric(x)))
  P.values <- unlist(lapply(Test_results, function(x) x$p.value))
  P.values <- as.data.frame(P.values)
  normality <- P.values > threshold  ##pavalues larger than 0.05 assume normal distribution 
  results <- cbind(P.values, normality)
  colnames(results) <- c("p.value", "normality")
  return(results)
}

##Levenes testing for variances
## dataframe is transposed with mz values as columns

#' @title levene test 
#' @description levene test for variance check
#' @param dataframe transposed dataframe with m/z as columns 
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold . the default is 0.05
#' @importFrom rstatix levene_test
leveneTestvalues <- function(dataframe, classifiers, threshold=0.05){
  Test_results <- apply(dataframe, 2, function(x) levene_test(formula = x ~ classifiers, data = dataframe))
  P.values <- unlist(lapply(Test_results, function(x) x$p))
  P.values <- as.data.frame(na.omit(P.values))
  Unequalvariances <- P.values < threshold
  results <- cbind(P.values, Unequalvariances)
  colnames(results) <- c("p.value", "unequal_variance")
  return(results)
}

#Welch T test
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers

#' @title welch t test 
#' @description welch t test between two groups of samples (assuming unequal variance)
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @param corr_method character string for correction method . default is "fdr"
#' @importFrom  future.apply future_apply
#' @importFrom miscTools colMedians
#' @importFrom utils combn
Welch_ttest <- function(dataframe, dataframe2, samples, threshold=0.1,corr_method="fdr"){
  mz <- as.numeric(row.names(t(dataframe)))
  T.test <- future_apply(dataframe, 2, function(x) t.test(x ~ samples, var.equal = FALSE))
  p.value <- unlist((lapply(T.test, function(x) x$p.value)))
  p.adj <- p.adjust(p.value, method = corr_method) 
  p.adj <- as.data.frame(p.adj)
  significant <- p.adj < threshold
  p.results <- cbind(p.value, p.adj, significant)
  colnames(p.results) <- c("p.value", "p.adj", "significant")
  Splitted_per_sampletype <- split(dataframe2, samples)
  Median_per_group <- lapply(Splitted_per_sampletype, function (x) colMedians(as.matrix(x)))
  Median_per_group <- as.data.frame(do.call(cbind, Median_per_group))
  combs <- combn(colnames(Median_per_group), 2)
  foldchange <- function(a, b) b/a
  foldchanges <- apply(combs, 2, function(col_names) foldchange(Median_per_group[, col_names[2]], Median_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = '_')
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  results <- cbind(mz, p.results, Median_per_group, foldchanges, log2foldchanges)
  return(results)
}


## wilcox test for non-parametric data
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers


#' @title wilcox test 
#' @description wilcox test between two groups of samples 
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @param corr_method character string for correction method . default is "fdr"
#' @param paired logical vector for paired test .. default paired = FALSE
#' @importFrom  future.apply future_apply
#' @importFrom miscTools colMedians
#' @importFrom utils combn
wilcox_results <- function(dataframe,dataframe2, samples, threshold=0.1,corr_method="fdr",paired = F){
  mz <- as.numeric(row.names(t(dataframe)))
  wilcox.test <- future_apply(dataframe, 2, function(x) wilcox.test(formula = x ~ samples,data=dataframe , paired = paired ))
  p.value <- unlist((lapply(wilcox.test, function(x) x$p.value)))
  p.adj <- p.adjust(p.value, method = corr_method) 
  p.adj <- as.data.frame(p.adj)
  significant <- p.adj < threshold
  p.results <- cbind(p.value, p.adj, significant)
  colnames(p.results) <- c("p.value", "p.adj", "significant")
  Splitted_per_sampletype <- split(dataframe2, samples)
  Median_per_group <- lapply(Splitted_per_sampletype, function (x) colMedians(as.matrix(x)))
  Median_per_group <- as.data.frame(do.call(cbind, Median_per_group))
  combs <- combn(colnames(Median_per_group), 2)
  foldchange <- function(a, b) b/a
  foldchanges <- apply(combs, 2, function(col_names) foldchange(Median_per_group[, col_names[2]], Median_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = '_')
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  results <- cbind(mz, p.results, Median_per_group, foldchanges, log2foldchanges)
  return(results)
}

## WelchAnova
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers



#' @title welch anova test 
#' @description welch anova test between more than two groups of samples (assuming unequal variance)
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @importFrom  future.apply future_apply
#' @importFrom miscTools colMedians
#' @importFrom utils combn
WelchAnova <- function(dataframe,dataframe2, samples, threshold=0.1){
  mz <- as.numeric(row.names(t(dataframe)))
  Test_results <- future_apply(dataframe, 2, function(x) oneway.test(x ~ samples, var.equal = FALSE))
  P.values <- unlist((lapply(Test_results, function(x) x$p.value)), use.names = FALSE)
  P.values <- as.data.frame(P.values)
  Significant_P.value <- P.values < threshold
  p.results <- cbind(P.values, Significant_P.value)
  colnames(p.results) <- c("p.value", "significant p.value")
  Splitted_per_sampletype <- split(dataframe2, samples)
  Mean_per_group <- lapply(Splitted_per_sampletype, function(x) colMeans(x))
  Mean_per_group <- as.data.frame(do.call(cbind, Mean_per_group))
  combs <- combn(colnames(Mean_per_group), 2)
  foldchange <- function(a, b) b/a
  foldchanges <- apply(combs, 2, function(col_names) foldchange(Mean_per_group[, col_names[2]], Mean_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = '_')
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  Median_per_group <- lapply(Splitted_per_sampletype, function (x) colMedians(x))
  Median_per_group <- as.data.frame(do.call(cbind, Median_per_group))
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  results <- cbind(mz, p.results, Median_per_group, foldchanges, log2foldchanges)
  return(results)
}




##Kruskall function
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers


#' @title Kruskal-Wallis test 
#' @description Kruskal-Wallis test  test between more than two groups of samples 
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @importFrom  future.apply future_apply
#' @importFrom miscTools colMedians
#' @importFrom utils combn
#' @importFrom dplyr filter
#' @importFrom dplyr filter_all
#' @importFrom dplyr any_vars
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
Kruskal_results <- function(dataframe,dataframe2,samples, threshold=0.1){
  mz <- (row.names(t(dataframe)))
  Test_results <- future_apply(dataframe, 2, function(x) kruskal.test(x, samples))
  p.values <- unlist((lapply(Test_results, function(x) x$p.value)), use.names = FALSE)
  p.values <- as.data.frame(p.values)
  Significant_p.value <- p.values < threshold
  p.results <- cbind(p.values, Significant_p.value)
  colnames(p.results) <- c("p.value", "significant p.value")
  Splitted_per_sampletype <- split(dataframe2, samples)
  Mean_per_group <- lapply(Splitted_per_sampletype, function(x) colMeans(x))
  Mean_per_group <- as.data.frame(do.call(cbind, Mean_per_group))
  combs <- combn(colnames(Mean_per_group), 2)
  foldchange <- function(a, b) b/a
  foldchanges <- apply(combs, 2, function(col_names) foldchange(Mean_per_group[, col_names[2]], Mean_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = '_')
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  Median_per_group <- lapply(Splitted_per_sampletype, function (x) colMedians(x))
  Median_per_group <- as.data.frame(do.call(cbind, Median_per_group))
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  results <- cbind(mz, p.results, Median_per_group, foldchanges, log2foldchanges)
  dataframe_t<-as.data.frame(t(dataframe))
  dataframe_t<-rownames_to_column(dataframe_t,"mz")
  sig_peaks_kruskall<-results %>% filter(results$`significant p.value`==TRUE)
  results_df<- dataframe_t %>% filter_all(any_vars(. %in% sig_peaks_kruskall$mz))
  results_df<- column_to_rownames(results_df,"mz")
  results_df<-as.data.frame(t(results_df))
  return(list("kruskal_results"=results,"significant_peaks"=results_df))
}

##Dunns test function
## dataframe is transposed dataframe with mz values as columns that is obtained from
## significant peaks from kruskall test (post hoc test )


#' @title Dunn test 
#' @description Dunn post hoc test after kruskal-wallis test
#' @param dataframe transposed dataframe with m/z as columns (results of significant kruskal test)
#' @param classifiers a factor vector with the classes of the samples 
#' @param threshold numerical value for wanted threshold 
#' @param comb a numerical value for number of group combinations
#' @param method a character string for correction method .. default is "bh"
#' @importFrom  future.apply future_apply
#' @importFrom dunn.test dunn.test
Dunn_test <- function(dataframe, samples, threshold,comb,method="bh"){
  mz <- rep(row.names(t(dataframe)), each = comb)
  Test_results <- future_apply(dataframe, 2, function(x) dunn.test(x, samples, method = method)) ## adjusted p value using bh
  P.adjusted <- as.data.frame(unlist((lapply(Test_results, function(x) x$P.adjusted)), use.names = FALSE))
  Significant_P.adjusted <- P.adjusted < threshold
  Groups <- as.data.frame(unlist((lapply(Test_results, function(x) x$comparisons)), use.names = FALSE))
  results <- cbind(mz, Groups, P.adjusted, Significant_P.adjusted)
  colnames(results) <-c("mz","Groups","P.adjusted", "Significant")
  return(results)
}


#' @title anova-tukey test
#' @description anova test followed by tukey HSD test as post-hoc
#' @param df transposed dataframe with m/z as columns 
#' @param classifiers a factor vector with the classes of the samples 
#' @importFrom  purrr map
anova_tukey<-function(df,classifiers){
  results <- map_dfr(df,
                   ~ TukeyHSD(aov( . ~ classifiers, data = df))[[1]][ ,4], .id = "mz")
  return(results)
  
  
}

