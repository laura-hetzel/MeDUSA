#' @title Shapiro's test for normality across cells
#' @description shapiro test for normality check
#' @param dataframe transposed dataframe with m/z as columns
#' @param assay
#' @importFrom stats shapiro.test
#' @export
shapiroTest <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(exp)
  dataframe <- assay(exp, assay)
  pvalues <- apply(dataframe, 1, function(x) shapiro.test(x)$p.value)
  rowData(exp)$shapiroTest <- pvalues
  exp
}

#' @title Levene test for unequal variances across phenotypes
#' @description levene test for variance check
#' @param exp transposed dataframe with m/z as columns
#' @param assay
#' @param classifiers a factor vector with the classes of the samples
#' @importFrom rstatix levene_test
#' @export
leveneTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype) {
  if (!validateExperiment(exp)) return(exp)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")

  data <- assay(exp, assay)
  classifiers <- as.factor(exp[[classifiers]])

  if (length(unique(classifiers)) > 2) {
    warning("Cannot perform Levene Test for phenotypes with > 2 groups")
    return(exp)
  }

  rowData(exp)$leveneTest <- apply(data, 1, function(x){
    levene_test(formula = x ~ classifiers, data = data)$p
  })
  exp
}

#' @title foldchange for two groups test
#' @param exp
#' @param assay
#' @param classifiers a factor vector with the classes of the samples
#' @importFrom utils combn
#' @export
foldChange <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype) {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")
  classifiers <- as.factor(exp[[classifiers]])
  if (length(unique(classifiers)) > 2) {
    warning("Cannot calculate Fold Change for phenotypes with > 2 groups")
    return(exp)
  }

  data <- as.data.frame(t(assay(exp, assay)))

  medians <- do.call(cbind, lapply(split(data, classifiers), colMedians))

  combs <- combn(colnames(medians), 2)

  foldchanges <- apply(combs, 2, function(names) {
    medians[, col_names[1]] / medians[, col_names[2]]
  })
  dimnames(foldchanges)[[2]] <- apply(combs, 2, paste, collapse = ".")
  rowData(exp)$foldChange <- as.data.frame(foldchanges)
  exp
}


# Welch T test
## dataframe is transposed with mz values as columns
## dataframe : scaled data
## dataframe2: raw intensities - imputed
## samples : classifiers

#' @title welch t test
#' @description welch t test between two groups of samples (assuming unequal variance)
#' @param exp
#' @param assay
#' @param classifiers
#' @param threshold numerical value for wanted threshold . the default is 0.1
#' @param corr_method character string for correction method . default is "fdr"
#' @importFrom stats p.adjust t.test
#' @export
welchTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                      threshold = 0.1, corr_method = "fdr") {
  if (!validateExperiment(exp)) return(NULL)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")

  classifiers <- as.factor(exp[[classifiers]])
  if (length(unique(classifiers)) > 2) {
    warning("Cannot perform t-test for phenotypes with > 2 groups")
    return(exp)
  }

  dataframe <- as.data.frame(t(assay(exp, assay)))

  p.value <- apply(dataframe, 2, function(x) t.test(x ~ classifiers, var.equal = FALSE)$p.value)
  p.adj <- p.adjust(p.value, method = corr_method)
  p.adj <- as.data.frame(p.adj)
  significant <- p.adj < threshold
  p.results <- cbind(p.value, p.adj, significant)
  colnames(p.results) <- c("p.value", "p.adj", "significant")
  rowData(exp)$welchTest <- p.results
  exp
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

  if (length(unique(classifiers)) > 2) {
    warning("Cannot perform Wilcoxon Test for phenotypes with > 2 groups")
    return(exp)
  }

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

#' @title welch anova test
#' @description welch anova test between more than two groups of samples (assuming unequal variance)
#' @param dataframe transposed dataframe with m/z as columns (transformed dataframe)
#' @param dataframe2 transposed dataframe with m/z as columns (raw imputed dataframe)
#' @param samples a factor vector with the classes of the samples
#' @param threshold numerical value for wanted threshold . the default is 0.1
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
#' @importFrom stats TukeyHSD aov
anova_tukey <- function(df, classifiers) {
  results <- do.call(rbind, lapply(seq_len(nrow(df)), function(i){
    TukeyHSD(aov(. ~ classifiers, data = df[i, ]))[[1]][, 4]
  }))
  return(results)
}
