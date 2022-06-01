




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
