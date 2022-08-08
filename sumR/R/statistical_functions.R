## shapiro's test for normality check
## dataframe is transposed with mz values as columns

#' @title shapiro test
#' @description shapiro test for normality check
#' @param dataframe transposed dataframe with m/z as columns
#' @param threshold numerical value for wanted threshold . the default is 0.05
#' @importFrom stats shapiro.test
#' @export
shapiroTest <- function(exp, assay = 1, threshold = 0.05) {
  if (!validateExperiment(exp)) return(NULL)
  dataframe <- as.data.frame(t(assay(exp, assay)))
  Test_results <- apply(dataframe, 2, function(x) shapiro.test(as.numeric(x)))
  P.values <- unlist(lapply(Test_results, function(x) x$p.value))
  P.values <- as.data.frame(P.values)
  normality <- P.values > threshold ## pavalues larger than 0.05 assume normal distribution
  results <- cbind(P.values, normality)
  colnames(results) <- c("p.value", "normality")
  rowData(exp)$shapiroTest <- results
  exp
}

## Levenes testing for variances
## dataframe is transposed with mz values as columns

#' @title levene test
#' @description levene test for variance check
#' @param dataframe transposed dataframe with m/z as columns
#' @param classifiers a factor vector with the classes of the samples
#' @param threshold numerical value for wanted threshold . the default is 0.05
#' @importFrom rstatix levene_test
#' @export
leveneTest <- function(exp, classifiers = metadata(exp)$phenotype,
                       assay = 1, threshold = 0.05, filter = FALSE) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))

  classifiers <- as.factor(exp[[classifiers]])
  Test_results <- apply(dataframe, 2, function(x) levene_test(formula = x ~ classifiers, data = dataframe))
  P.values <- unlist(lapply(Test_results, function(x) x$p))
  P.values <- as.data.frame(na.omit(P.values))
  Unequalvariances <- P.values < threshold
  results <- cbind(P.values, Unequalvariances)
  colnames(results) <- c("p.value", "unequal_variance")
  rowData(exp)$leveneTest <- results
  if (any(results$unequal_variance)) {
    if (filter) exp <- exp[results$unequal_variance, ]
  }
  exp
}

#' @title foldchange for two groups test
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param samples a factor vector with the classes of the samples
#' @importFrom miscTools colMedians
#' @importFrom utils combn
#' @export
foldChange <- function(exp, classifiers = metadata(exp)$phenotype, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe2 <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])
  Splitted_per_sampletype <- split(dataframe2, classifiers)
  Median_per_group <- do.call(cbind, lapply(Splitted_per_sampletype, function(x) colMedians(as.matrix(x))))
  combs <- combn(colnames(Median_per_group), 2)
  fc <- function(a, b) b / a
  foldchanges <- apply(combs, 2, function(col_names) fc(Median_per_group[, col_names[2]], Median_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = "_")
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  results <- cbind(Median_per_group, foldchanges, log2foldchanges)
  rowData(exp)$foldChange <- results
  exp
}




#' @title foldchange for more than two groups test
#'
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param samples a factor vector with the classes of the samples
#' @importFrom miscTools colMedians
#' @importFrom utils combn
foldChangeGroups <- function(dataframe2, samples) {
  Splitted_per_sampletype <- split(dataframe2, samples)
  Mean_per_group <- do.call(cbind, lapply(Splitted_per_sampletype, function(x) colMeans(x)))
  combs <- combn(colnames(Mean_per_group), 2)
  fc <- function(a, b) b / a
  foldchanges <- apply(combs, 2, function(col_names) fc(Mean_per_group[, col_names[2]], Mean_per_group[, col_names[1]]))
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = "_")
  foldchanges <- as.data.frame(foldchanges)
  colnames(foldchanges) <- paste("fc", colnames(foldchanges), sep = "_")
  log2foldchanges <- log2(foldchanges)
  colnames(log2foldchanges) <- paste("log2", colnames(foldchanges), sep = "")
  Median_per_group <- lapply(Splitted_per_sampletype, function(x) colMedians(x))
  Median_per_group <- as.data.frame(do.call(cbind, Median_per_group))
  colnames(Median_per_group) <- paste("median", colnames(Median_per_group), sep = "_")
  results <- cbind(Median_per_group, foldchanges, log2foldchanges)
  return(results)
}



# Welch T test
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers

#' @title welch t test
#' @description welch t test between two groups of samples (assuming unequal variance)
#'
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @param samples a factor vector with the classes of the samples
#' @param corr_method character string for correction method . default is "fdr"
#' @importFrom stats p.adjust t.test
#' @export
welchTest <- function(exp, classifiers = metadata(exp)$phenotype, assay = 1, threshold = 0.1,
                      corr_method = "fdr") {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])
  p.value <- apply(dataframe, 2, function(x) t.test(x ~ classifiers, var.equal = FALSE)$p.value)
  p.adj <- p.adjust(p.value, method = corr_method)
  p.adj <- as.data.frame(p.adj)
  significant <- p.adj < threshold
  p.results <- cbind(p.value, p.adj, significant)
  colnames(p.results) <- c("p.value", "p.adj", "significant")
  rowData(exp)$welchTest <- p.results
  exp
}

#' @title Keep the most variable features
#' @param exp
#' @param assay
#' @param top
#' @export
keepVariableFeatures <- function(exp){
  if (!validateExperiment(exp)) return(NULL)
  if (!"leveneTest" %in% colnames(rowData(exp))) {
    message("'leveneTest' not found in rowData. Please run 'leveneTest' first.")
    return(exp)
  }
  exp[which(rowData(exp)$leveneTest$unequal_variance), ]
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
#' @param samples a factor vector with the classes of the samples
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @param corr_method character string for correction method . default is "fdr"
#' @param paired logical vector for paired test .. default paired = FALSE
#' @importFrom stats p.adjust wilcox.test
wilcoxTest <- function(exp, classifiers = metadata(exp)$phenotype,
                       assay = 1, threshold = 0.1,
                       corr_method = "fdr", paired = F) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])

  wilcox.test <- apply(dataframe, 2, function(x) wilcox.test(formula = x ~ classifiers,
                                                             data = dataframe, paired = paired))
  p.value <- unlist((lapply(wilcox.test, function(x) x$p.value)))
  p.adj <- p.adjust(p.value, method = corr_method)
  p.adj <- as.data.frame(p.adj)
  significant <- p.adj < threshold
  p.results <- cbind(p.value, p.adj, significant)
  colnames(p.results) <- c("p.value", "p.adj", "significant")
  rowData(exp)$wilcoxTest <- p.results
  exp
}

## WelchAnova
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers



#' @title welch anova test
#' @description welch anova test between more than two groups of samples (assuming unequal variance)
#'
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param samples a factor vector with the classes of the samples
#' @param threshold numerical value for wanted threshold . the default is 0.1
#'
#' @importFrom stats oneway.test
welchAnova <- function(exp, classifiers = metadata(exp)$phenotype,
                       assay = 1, threshold = 0.1) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])

  P.values <- apply(dataframe, 2, function(x) oneway.test(x ~ classifiers, var.equal = FALSE)$p.value)
  Significant_P.value <- P.values < threshold
  p.results <- cbind(P.values, Significant_P.value)
  colnames(p.results) <- c("p.value", "significant p.value")
  rowData(exp)$welchAnova <- p.results
  exp
}




## Kruskall function
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers


#' @title Kruskal-Wallis test
#' @description Kruskal-Wallis test  test between more than two groups of samples
#'
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param samples a factor vector with the classes of the samples
#' @param threshold numerical value for wanted threshold . the default is 0.1
#'
#' @importFrom stats kruskal.test
#' @export
kruskalTest <- function(exp, classifiers = metadata(exp)$phenotype,
                        assay = 1, threshold = 0.1) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])
  p.values <- apply(dataframe, 2, function(x) kruskal.test(x, classifiers)$p.value)
  Significant_p.value <- p.values < threshold
  p.results <- cbind(p.values, Significant_p.value)
  colnames(p.results) <- c("p.value", "significant p.value")
  rowData(exp)$kruskalTest <- p.results
  exp
}

## Dunns test function
## dataframe is transposed dataframe with mz values as columns that is obtained from
## significant peaks from kruskall test (post hoc test )


#' @title Dunn test
#' @description Dunn post hoc test after kruskal-wallis test
#'
#' @param dataframe transposed dataframe with m/z as columns (results of significant kruskal test)
#' @param threshold numerical value for wanted threshold default is 0.05
#' @param comb a numerical value for number of group combinations (4 groups = 6 combination)
#' @param samples a factor vector with the classes of the samples
#' @param method a character string for correction method .. default is "bh"
#'
#' @importFrom dunn.test dunn.test
#' @export
dunnTest <- function(exp, classifiers = metadata(exp)$phenotype,
                     assay = 1, threshold = 0.05, comb = 6, method = "bh") {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- as.factor(exp[[classifiers]])

  Test_results <- apply(dataframe, 2, function(x) dunn.test(x, classifiers, method = method)) ## adjusted p value using bh
  P.adjusted <- as.data.frame(unlist((lapply(Test_results, function(x) x$P.adjusted)), use.names = FALSE))
  Significant_P.adjusted <- P.adjusted < threshold
  Groups <- as.data.frame(unlist((lapply(Test_results, function(x) x$comparisons)), use.names = FALSE))
  results <- cbind(Groups, P.adjusted, Significant_P.adjusted)
  colnames(results) <- c("Groups", "P.adjusted", "Significant")
  rowData(exp)$dunnTest <- results
  exp
}


#' @title anova-tukey test
#' @description anova test followed by tukey HSD test as post-hoc
#' @param df transposed dataframe with m/z as columns
#' @param classifiers a factor vector with the classes of the samples
#' @importFrom  purrr map
anova_tukey <- function(df, classifiers) {
  results <- map_dfr(df,
    ~ TukeyHSD(aov(. ~ classifiers, data = df))[[1]][, 4],
    .id = "mz"
  )
  return(results)
}
