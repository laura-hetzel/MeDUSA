#' @title Shapiro's test for normality across cells
#' @description shapiro test for normality check
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param method
#' @importFrom stats shapiro.test
#' @export
shapiroTest <- function(exp, assay = 1, method = "fdr") {
  if (!validateExperiment(exp)) return(exp)
  dataframe <- assay(exp, assay)
  p.values <- apply(dataframe, 1, function(x) shapiro.test(x)$p.value)
  rowData(exp)$shapiroTest <- p.adjust(p.values, method = method)
  exp
}

#' @title foldchange for two groups test
#' @description
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers a factor vector with the classes of the samples
#' @importFrom utils combn
#' @export
foldChange <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype) {
  if (!testCheck(exp, classifiers, maxGroups = 2)) return(exp)

  classifiers <- exp[[classifiers]]

  first <- exp[, classifiers == unique(classifiers)[1]]
  meds1 <- rowMedians(as.matrix(assay(first, 1)))

  second <- exp[, classifiers == unique(classifiers)[2]]
  meds2 <- rowMedians(as.matrix(assay(second, 1)))

  rowData(exp)$foldChange <- meds1 / meds2
  exp
}

#' @title Levene test for unequal variances across phenotypes
#' @description levene test for variance check
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @param method
#' @importFrom rstatix levene_test
#' @export
leveneTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                       method = "fdr") {
  if (!testCheck(exp, classifiers, maxGroups = 2)) return(exp)

  classifiers <- as.factor(exp[[classifiers]])
  data <- as.data.frame(t(assay(exp, assay)))
  p.values <- apply(data, 2, function(x){
    rstatix::levene_test(formula = x ~ classifiers, data = data)$p
  })
  rowData(exp)$leveneTest <- p.adjust(p.values, method = method)
  exp
}


#' @title welch t test
#' @description welch t test between two groups of samples (assuming unequal
#' variance)
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @param method
#' @importFrom stats p.adjust t.test
#' @export
welchTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                      method = "fdr") {
  if (!testCheck(exp, classifiers, maxGroups = 2)) return(exp)

  classifiers <- as.factor(exp[[classifiers]])
  dataframe <- as.data.frame(t(assay(exp, assay)))
  p.value <- apply(dataframe, 2, function(x){
    t.test(x ~ classifiers, var.equal = FALSE)$p.value
  })
  rowData(exp)$welchTest <- p.adjust(p.value, method = method)
  exp
}

#' @title wilcox test
#' @description wilcox test between two groups of samples
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @param method
#' @importFrom stats p.adjust wilcox.test
#' @export
wilcoxTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                       method = "fdr") {

  if (!testCheck(exp, classifiers, maxGroups = 2)) return(exp)

  classifiers <- as.factor(exp[[classifiers]])
  dataframe <- as.data.frame(t(assay(exp, assay)))
  wilcox.test <- apply(dataframe, 2, function(x){
    wilcox.test(formula = x ~ classifiers, data = dataframe, paired = FALSE)
  })
  p.value <- unlist((lapply(wilcox.test, function(x) x$p.value)))
  rowData(exp)$wilcoxTest <- p.adjust(p.value, method = method)
  exp
}

#' @title welch anova test
#' @description welch anova test between more than two groups of samples
#' (assuming unequal variance)
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @param method
#' @importFrom stats oneway.test
#' @export
welchAnova <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                       method = "fdr") {
  if (!testCheck(exp, classifiers, minGroups = 3)) return(exp)

  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- exp[[classifiers]]

  p.values <- apply(dataframe, 2, function(x){
    oneway.test(formula = x ~ classifiers, data = dataframe,
                subset = x)$p.value
  })
  rowData(exp)$welchAnova <- p.adjust(p.values, method = method)
  exp
}

#' @title Dunn test
#' @description Dunn post hoc test after kruskal-wallis test
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @param method a character string for correction method .. default is "bh"
#' @importFrom dunn.test dunn.test
#' @export
dunnTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                     method = "bh") {

  if (!testCheck(exp, classifiers, minGroups = 3)) return(exp)

  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- exp[[classifiers]]

  p.value <- lapply(seq_len(ncol(dataframe)),function(i){
    test <- quiet(dunn.test(dataframe[,i], classifiers, method = method))
    setNames(test$P.adjust, test$comparisons)
  })

  rowData(exp)$dunnTest <- as.data.frame(do.call(rbind, p.value))
  exp
}

#' @title anova-tukey test
#' @description anova test followed by tukey HSD test as post-hoc
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @importFrom stats TukeyHSD aov
#' @export
aovTukeyTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype) {
  if (!testCheck(exp, classifiers, minGroups = 3)) return(exp)

  dataframe <- as.data.frame(t(assay(exp, assay)))
  classifiers <- exp[[classifiers]]

  res <- lapply(seq_len(ncol(dataframe)), function(i){
    df <- data.frame(points = dataframe[, i], classifiers = classifiers)
    res <- TukeyHSD(aov(formula("points ~ classifiers"), data = df))
    res$classifiers[, "p adj"]
  })

  rowData(exp)$aovTukeyTest <- as.data.frame(do.call(rbind, res))
  exp
}

#' @title Check if the test picked can be executed
#' @description
#' @details
#' @returns
#' @param exp
#' @param classifiers
#' @param minGroups
#' @param maxGroups
testCheck <- function(exp, classifiers, minGroups = 0, maxGroups = Inf) {
  if (!validateExperiment(exp)) return(FALSE)
  if (is.null(classifiers)) stop("Cannot perform test without classifiers")

  groups <- length(unique(exp[[classifiers]]))
  if (groups >= minGroups & groups <= maxGroups) {
    return(TRUE)
  }
  warning("Inadequate number of phenotypes in the given phenotype column")
  return(FALSE)
}

#' @title Detect an appropriate statistical test
#' @description
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @export
detectTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                       threshold = 2/3){
  groups <- length(unique(exp[[classifiers]]))
  exp <- shapiroTest(exp, assay)
  type <- ifelse(sum(rowData(exp)$shapiroTest > 0.05) / nrow(exp) > threshold,
         "normal", "notNormal")
  if (groups < 3 & type == "normal") {
    return("welchTest")
  } else if (groups < 3 & type == "notNormal") {
    return("wilcoxTest")
  } else if (groups >= 3 & type == "normal") {
    return("aovTukeyTest")
  }
  return("dunnTest")
}

#' @title Detect and execute an appropriate statistical test
#' @description
#' @details
#' @returns
#' @param exp
#' @param assay
#' @param classifiers
#' @export
autoTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype){
  test <- detectTest(exp, assay, classifiers)
  message("Using test: ", test)
  get(test)(exp, assay, classifiers)
}
