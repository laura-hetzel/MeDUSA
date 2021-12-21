library(plyr)
library(future.apply) #statistical functions 
library(robustbase) ## col mean and median in statisitical tests 
library(tibble)
library(dunn.test) ## dunn test 
library(data.table)
library(miscTools)
library(rstatix)


# statistical tests  ----------------------------------------------------

##Levenes testing for variances
## dataframe is transposed with mz values as columns
leveneTestvalues <- function(dataframe, samples, threshold){
  Test_results <- future_apply(dataframe, 2, function(x) levene_test(formula = x ~ classifiers, data = dataframe))
  P.values <- unlist(lapply(Test_results, function(x) x$p))
  P.values <- as.data.frame(na.omit(P.values))
  Unequalvariances <- P.values < threshold
  results <- cbind(P.values, Unequalvariances)
  colnames(results) <- c("p.value", "unequalvariance")
  return(results)
}

## shapiro's test for normality check
## dataframe is transposed with mz values as columns

shapiro_results <- function(dataframe,threshold){
  Test_results <- future_apply(dataframe,2, function(x) shapiro.test(as.numeric(x)))
  P.values <- unlist(lapply(Test_results, function(x) x$p.value))
  P.values <- as.data.frame(P.values)
  normality <- P.values > threshold  ##pavalues larger than 0.05 assume normal distribution 
  results <- cbind(P.values, normality)
  colnames(results) <- c("p.value", "normality")
  return(results)
}



##Welch T test
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers
Welchttest <- function(dataframe, dataframe2, samples, threshold){
  mz <- as.numeric(row.names(t(dataframe)))
  T.test <- future_apply(dataframe, 2, function(x) t.test(x ~ samples, var.equal = FALSE))
  p.value <- unlist((lapply(T.test, function(x) x$p.value)))
  p.adj <- p.adjust(p.value, method = "fdr") 
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
wilcox_results <- function(dataframe,dataframe2, samples, threshold){
  mz <- as.numeric(row.names(t(dataframe)))
  wilcox.test <- future_apply(dataframe, 2, function(x) wilcox.test(formula = x ~ samples,data=dataframe , paired = F ))
  p.value <- unlist((lapply(wilcox.test, function(x) x$p.value)))
  p.adj <- p.adjust(p.value, method = "fdr") 
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
WelchAnova <- function(dataframe,dataframe2, samples, threshold){
  mz <- as.numeric(row.names(t(dataframe)))
  Test_results <- future_apply(dataframe, 2, function(x) oneway.test(x ~ samples, var.equal = FALSE))
  P.values <- unlist((lapply(Test_results, function(x) x$p.value)), use.names = FALSE)
  P.values <- as.data.frame(P.values)
  Significant_P.value <- P.values < threshold
  p.results <- cbind(P.values, Significant_P.value)
  colnames(p.results) <- c("p.value", "significant p.value?")
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
Kruskall <- function(dataframe,dataframe2,samples, threshold){
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
  return(results)
}

##Dunns test function
## dataframe is transposed dataframe with mz values as columns that is obtained from
## significant peaks from kruskall test (post hoc test )
Dunn_test <- function(dataframe, samples, threshold){
  mz <- rep(row.names(t(dataframe)), each = 6) ## number = 6 depends on the number of groups and combination
  # here we have 4 groups so 6 combination
  Test_results <- future_apply(dataframe, 2, function(x) dunn.test(x, classifiers, method = "bh")) ## adjusted p value using bh
  P.adjusted <- as.data.frame(unlist((lapply(Test_results, function(x) x$P.adjusted)), use.names = FALSE))
  Significant_P.adjusted <- P.adjusted < threshold
  Groups <- as.data.frame(unlist((lapply(Test_results, function(x) x$comparisons)), use.names = FALSE))
  results <- cbind(mz, Groups, P.adjusted, Significant_P.adjusted)
  colnames(results) <-c("mz","Groups","P.adjusted", "Significant?")
  return(results)
}


## Anova followed by tukey hsd : post hoc test after anova for parametric data (more than two groups)
## dataframe is transposed with mz values as columns
anova <- map_dfr(master_df_t,
                 ~ TukeyHSD(aov( . ~ classifiers, data = master_df))[[1]][ ,4], .id = "mz")
