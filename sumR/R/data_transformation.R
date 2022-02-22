
#' @title data imputation
#' @description This function apply data imputation from 1 to selected noise level  
#' @param data dataframe 
#' @param noise numerical value for the noise level
#' @param seed global seed for reproducible results
#' @importFrom stats runif
  data_imputation<-function(data,noise,seed){
    set.seed(seed)
    data[data == 0] <- runif(1, min = 1, max = noise)
    return(data)
}
  



#' @title data transformation
#' @description This function apply data transformation using log2 and pareto scale 
#'
#' @param centering 
#' @param data dataframe 
#'
#' @importFrom tibble column_to_rownames
#' @importFrom IMIFA pareto_scale
data_transform<-function(data,centering=F){
  
  data<-column_to_rownames(data,"mz")
  data<-log2(data)
  data<-pareto_scale((data),centering = centering)
  data<-as.data.frame(data)
  return(data)
  
}
