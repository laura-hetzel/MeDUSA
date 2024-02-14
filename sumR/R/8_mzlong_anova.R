
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
mzlong_analysis_anova <- function(input_mzlong_obj, metadata, phenotypes, p_cutoff = 0., cores = 2,){
  metadata <- local.meta_polarity_fixer(input_mzlong_obj, metadata)  
  mzlong <- dplyr::left_join(input_mzlong_obj, metadata[c("sample_name", "phenotype")], by = "sample_name")
  mzlong<-mzlong[mzlong$phenotype %in% phenotypes,]
  uniq_mz <- unique(mzlong$mz)
  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    #TODO: should mz be considered a factor?
    #anova <- aov(log ~ mz + phenotype , data = mzlong)
    anova <- pbapply::pblapply(uniq_mz, anova.by_mz, mzlong=mzlong)
    summary <- pbapply::pblapply(seq_along(anova), anova.summary, anova=anova, mzs = uniq_mz)
    do.call(rubind, summary)
    tukey <- pbapply::pblapply(seq_along(anova), anova.get_p, anova = anova, mzs = uniq_mz)
    tukey <- do.call(rbind,tukey)
  }, finally={
    local.kill_threads(cl)
  })
  summary <- summary[order(-summary$`Pr(>F)`, decreasing = T),]
  imp_mz <- summary$mz[1:100]
  #best <- AICcmodavg::aictab(a,modnames = uniq_mz)
  #best_index <- match(best$Modnames, uniq_mz)
  list("anova" = anova, "summary" = summary, "tukey" = tukey, "imp_mz" = imp_mz)
}

anova.by_mz <- function(mz, mzlong){
  input_mz <- mzlong[ mzlong$mz == mz,]
  #TODO make sure + not * (are the vars independant or not?)
  aov(log ~ phenotype , data = input_mz)
}

anova.get_p <- function(index, anova, mzs){
  tmp <- as.data.frame(TukeyHSD(anova[[index]])$phenotype)
  tmp$mz <- mzs[index]
  tmp
}

anova.summary <- function(index, anova, mzs){
  tmp <- summary(anova[[index]])[[1]]["phenotype",]
  tmp$mz <- mzs[index]
  tmp
}