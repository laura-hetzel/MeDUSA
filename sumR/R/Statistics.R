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
#' @description This function calculates the median fold change between the
#' phenotypes determined by `classifiers`.
#' @details The fold change is a metric for the difference between groups and
#' can be indicative of the significance of the differences. When the assay has
#' been scaled with a log transformation, this effectively calculates the logFC,
#' often used for volcano plots.
#' @returns SummarizedExperiment with `foldChange` column in the [rowData] slot.
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
#' @description This function tests the variance across phenotypes for each
#' peak and adds the corresponding p-value to the column `leveneTest` to the
#' [rowData] slot.
#' @details The levene test is a test for calculating unequal variances across
#' groups. This can be an indicator for peaks that have different values in the
#' provided assay.
#' @returns SummarizedExperiment with `leveneTest` column in the [rowData] slot.
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


#' @title Welch T-test for 2 groups
#' @description Welch T-test between two groups of samples (assuming unequal
#' variance).
#' @details The Welch T-test is a two-sample test to see if the two
#' phenotypes have equal means in their distribution. This test assums unequal
#' variance across the phenotypes.
#' @returns SummarizedExperiment with `welchTest` column in the [rowData] slot.
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

#' @title Wilcoxon test for non-parametric distributed data of 2 groups
#' @description This function will perform the Wilcoxon signed-rank test for
#' non-parametric data. The groups are determined by the phenotype given and
#' will produce an error if the number of groups is not equal to 2.
#' @details This test is a common test that tests the difference between the
#' phenotypes proivded. It is considered the equivalent to the [welchTest], but
#' for non-parametric data, which causes it to be less powerful in general.
#' @returns SummarizedExperiment with `wilcoxTest` column in the [rowData] slot.
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
#' @returns SummarizedExperiment with `dunnTest` column in the [rowData] slot.
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

#' @title Anova with Tukey HSD post-hoc test for > 2 groups
#' @description This function excutes an anova test and consecutively the
#' Tukey HSD post-hoc test for determining which compounds are significantly
#' different across > 2 phenotypes.
#' @details
#' @returns SummarizedExperiment with `aovTukeyTest` column in the [rowData]
#' slot.
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
#' @description This function checks if the test is suitable for the data
#' that is provided and will throw a warning if not or and error when no
#' classifiers are provided.
#' @details This function is called by all tests to see if the test is
#' appropriate as a check. Used internally by sumR.
#' @returns Boolean value, is the test valid to use?
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
#' @description This function checks the supplied assay and phenotypes for
#' normality and number of groups respectively. Consequently, it picks an
#' appropriate test to check for statistical difference.
#' @details This function uses the [shapiroTest] to determine for each peak
#' if the assay supplied is normally distributed. Together with the number of
#' groups found in the phenotype slot set by [setPhenotype], 4 different tests
#' can be suggested:
#'
#' - __welchTest__: Number of groups equals 2 and data is considered normally
#' distributed
#' - __wilcoxTest__: Number of groups equals 2 and data is considered not to be
#' normally distributed.
#' - __aovTukeyTest__: Number of groups is higher than 2 and the data is
#' considered normally distributed.
#' - __dunnTest__: Number of groups is higher than 2 and the data is
#' considered not to be normally distributed.
#'
#' @returns A character vector with the name of the test that is suggested.
#' This is equal to the function name.
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param threshold What fraction of the peaks need to be normally distributed
#' before performing a parametric statistical test? Defaults to `2/3`.
#' @export
#' @seealso [autoTest]
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
#' @description This function checks the supplied assay and phenotypes for
#' normality and number of groups respectively. Consequently, it picks an
#' appropriate test to check for statistical difference. This function will
#' print a message with the test chosen. See Details for more explanation on
#' which tests are considered.
#' @inherit detectTest details
#' @returns SummarizedExperiment with an added column in [rowData] with the
#' results of the chosen test.
#' @param exp SummarizedExperiment with a (auto)scaled assay.
#' @param assay Name or index of the assay with (auto)scaled values. Defaults
#' to the first assay (index 1).
#' @param classifiers Column in the [colData] that describes the phenotype or
#' predictor class used for statistical testing and modelling. Defaults to the
#' phenotype set using [setPhenotype].
#' @param threshold Used for [detectTest]. What fraction of the peaks need to
#' be normally distributed before performing a parametric statistical test?
#' Defaults to `2/3`.
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
autoTest <- function(exp, assay = 1, classifiers = metadata(exp)$phenotype,
                     threshold = 2/3){
  test <- detectTest(exp, assay, classifiers, threshold)
  message("Using test: ", test)
  get(test)(exp, assay, classifiers)
}
