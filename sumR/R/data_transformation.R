

#' @title data imputation
#' @description This function apply data imputation from 1 to selected noise level  
#' @param data dataframe 
#' @param noise numerical value for the noise level
#' @param seed global seed for reproducible results
#' @importFrom  purrr map_dfc
data_imputation<-function(data,noise,seed){
  
  map_dfc(data, function(x){
    set.seed(seed)
    for (i in 1:length(x)) {
      if(x[i] == 0) {
        x[i] <- runif(1, min = 1, max = noise)
      }
      else {next}
    }
    return(x)
  })
  
}


#' @title data transformation
#' @description This function apply data transformation using log2 and pareto scale 
#' @param data dataframe 
#' @param centring a logical vector TRUE or FALSE for centered scaling . the default is FALSE
#' @importFrom tibble column_to_rownames
#' @importFrom IMIFA pareto_scale
data_transform<-function(data,centering=F){
  
  data<-column_to_rownames(data,"mz")
  data<-log2(data)
  data<-pareto_scale((data),centering = centering)
  data<-as.data.frame(data)
  return(data)
  
}