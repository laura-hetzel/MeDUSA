
## choosing alpha value for elastic net
#' @title Alpha value choice
#' @description This function choose the alpha value for elastic net regression with the desired accuracy 
#' @param training training set
#' @param test test set
#' @param seed global seed for reproducible results
#' @param type.measure loss to use for cross-validation. type.measure="class" ,"auc"
#' @importFrom glmnet cv.glmnet
#' @importFrom caret confusionMatrix
#' @importFrom dplyr select
glmnet_cv_alphachoice <- function(training,test,type.measure,seed,lambda.min){
  
  list_fits<-list()
  for(i in 0:10){
    
    fit_name <-paste0("alpha",i/10)
    set.seed(seed)
    list_fits[[fit_name]] <-  cv.glmnet(y = training$samples,
                                        x = as.matrix(scale(training %>% dplyr::select(-samples))),
                                        family = "binomial",alpha=i/10,type.measure = type.measure)
    
  }
  
  results <-data.frame()
  for(i in 0:10){
    fit_name <-paste0("alpha",i/10)
    predict <-predict(list_fits[[fit_name]],s=list_fits[[fit_name]]$lambda.min,newx = as.matrix(scale(test %>% dplyr::select(-samples))),
                      type="class")
    cf<-confusionMatrix(table(data=predict, reference = test$samples))
    cf_acc<-cf[["overall"]][["Accuracy"]]
    temp <-data.frame(alpha=i/10,cf_acc=cf_acc,fit_name=fit_name)
    results<-rbind(results,temp)
    
    
  }
  return(results)
  
}


#
#' @title final glmnet model
#' @description This function makes the final model of glmnet using the desired alpha value 
#' @param training training set
#' @param test test set
#' @param seed global seed for reproducible results
#' @param type.measure loss to use for cross-validation. type.measure="class" ,"auc"
#' @param master_df the master df of the raw or scaled intensitites that have m/z column
#' @importFrom glmnet cv.glmnet
#' @importFrom caret confusionMatrix
#' @importFrom dplyr select
#' @importFrom dplyr filter_all
#' @importFrom dplyr any_vars
final_glmnet <- function(training,set,alpha,type.measure,master_df,seed){
  set.seed(seed)
  fit_final <- cv.glmnet(y = training$samples,
                         x = as.matrix(scale(training %>% dplyr::select(-samples))),
                         family = "binomial",alpha=alpha,type.measure = type.measure)
  plot(fit_final)
  print(fit_final)
  
  predict <-predict(fit_final,s="lambda.min",newx = as.matrix(scale((test %>% dplyr::select(-samples)))),type="class")
  print(confusionMatrix(table(data=predict, reference = test$samples)))
  
  k <-  coef(fit_final, s = "lambda.min"); k2 <-  k[which(k!=0),]; k2
  
  k<-as.data.frame(names(k2)) ; k<-k[-1,]
  imp_var<-master_df%>% filter_all(any_vars(. %in% k)) ## obtaining the intensities of the all imp mz from the master df
  
  
  return(imp_var)
}
