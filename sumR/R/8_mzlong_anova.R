# *** Anova Analysis -----------------------------------------------------
#' Anova my shiz
#'
#' @param input_mzlong_obj \cr
#'   DataFrame : LongLog2 of Input MZ-Obj
#' @param metadata \cr
#'   DataFrame :meta
#' @param attribute \cr
#'   String: Which column in metadata matters?
#'
#' @returns 
#'
#' @export
mzlong_analysis_anova <- function(input_mzlong_obj, metadata, attribute = "phenotype", cores = 2){
  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    input_mzlong_obj$att <- pbapply::pbapply(input_mzlong_obj,1,cl=cl,function(i){  
      metadata[metadata$sample_name == input_mzlong_obj[i,]["sample"]][attribute]
    })
    anova <- aov(mz ~ att, data = input_mzlong_obj)
  
    
  },
  finally={
    #This doesn't seem to stop the cluster :/
    local.kill_threads(cl)
  })
  return(anova)
}