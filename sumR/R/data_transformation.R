#' @title data imputation
#' @description This function apply data imputation from 1 to selected noise level
#' @param data dataframe
#' @param noise numerical value for the noise level
#' @param seed global seed for reproducible results
#' @importFrom stats runif
data_imputation <- function(data, noise, seed=42) {
  set.seed(seed)
  data[data == 0] <- runif(length(data == 0), min = 1, max = noise)
  return(data)
}




#' @title data transformation
#' @description This function apply data transformation using log2 and pareto scale
#' @param data dataframe
#' @importFrom tibble column_to_rownames
data_transform <- function(data) {
  data <- column_to_rownames(data, "mz")
  data <- log2(data)
  data <- paretoScaling((data))
  data <- as.data.frame(data)
  return(data)
}



#' @title pareto scaling
#' @param df dataframe
paretoScaling <- function(df) apply(df, 2, function(x) (x - mean(x)) / sqrt(sd(x)))
