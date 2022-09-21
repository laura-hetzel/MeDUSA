#' @title Shapiro's test for normality across cells
#' @description Shapiro's test for normality is a test that indicates if the
#' supplied distribution is significantly different from a normal distribution.
#' This test stored the results of the adjusted p-value in the [rowData] slot
#' of the `exp`.
#' @details A significant adjusted p-value indicates that the distribution is
#' very likely not normal, and therefore non-parametric tests should be used to
#' get an accurate indication of significance, e.g. [wilcoxTest] for 2 groups
#' and [dunnTest] for > 2 groups.
#'
#' It is recommended to impute an assay using e.g. [imputation] and scale it
#' using [autoScale] to get a fair estimation of normality.
#' @returns SummarizedExperiment with updated [rowData] slot. The results will
#' be stored in the `shapiroWilk` column of the rowData of `exp`.
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param method Correction method to address the multiple testing problem.
#' Defaults to `"fdr"` (False Discovery Rate).
#' @importFrom stats shapiro.test
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Shapiro's test
#' sumRnegative <- shapiroTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$shapiroTest
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Fold Change
#' sumRnegative <- foldChange(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$foldChange
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param method Correction method to address the multiple testing problem.
#' Defaults to `"fdr"` (False Discovery Rate).
#' @importFrom rstatix levene_test
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Levene's Test
#' sumRnegative <- leveneTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$leveneTest
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param method Correction method to address the multiple testing problem.
#' Defaults to `"fdr"` (False Discovery Rate).
#' @importFrom stats p.adjust t.test
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Welch' T-Test
#' sumRnegative <- welchTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$welchTest
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param method Correction method to address the multiple testing problem.
#' Defaults to `"fdr"` (False Discovery Rate).
#' @importFrom stats p.adjust wilcox.test
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Wilcoxon Test
#' sumRnegative <- wilcoxTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$wilcoxTest
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

#' @title Dunn test
#' @description Dunn post hoc test after kruskal-wallis test
#' @details
#' @returns
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param Correction method to address the multiple testing problem.
#' Defaults to `"bh"` (Benjamini-Hochberg).
#' @importFrom dunn.test dunn.test
#' @seealso [dunn.test] for correction methods available.
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Dunn's Test
#' sumRnegative <- dunnTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$dunnTest
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @importFrom stats TukeyHSD aov formula
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Run Anova & Tukey Test
#' sumRnegative <- aovTukeyTest(sumRnegative)
#'
#' # Show results of first 5 compounds
#' rowData(sumRnegative[1:5, ])$aovTukeyTest
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param minGroups Minimum number of groups for the test. Defaults to `0`
#' @param maxGroups Maximum number of groups for the test. Defaults to `Inf`
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Detect which test is suitable
#' detectTest(sumRnegative)
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
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @export
#' @examples
#' #' # Read example data
#' data("sumRnegative")
#'
#' # Set colData and phenotype
#' df <- read.csv(system.file("cellData.csv", package = "sumR"), row.names = 1)
#' sumRnegative <- addCellData(sumRnegative, df)
#' sumRnegative <- setPhenotype(sumRnegative, "Treatment")
#'
#' # Imputate and Scale data
#' sumRnegative <- imputation(sumRnegative)
#' sumRnegative <- autoScale(sumRnegative)
#'
#' # Detect and execute appropriate test
#' sumRnegative <- autoTest(sumRnegative)
autoTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype){
  test <- detectTest(exp, assay, classifiers)
  message("Using test: ", test)
  get(test)(exp, assay, classifiers)
}
