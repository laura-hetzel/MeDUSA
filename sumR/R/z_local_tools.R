### z_* written by ehetzel zz_* written by others

local.dir_sep <- function() {
  if ( .Platform$OS.type == "unix" ){
    dir_sep = "/" #Unix
  } else {
    dir_sep = "\\" #Windows, note escape char
  }
  dir_sep
}

#Show zeros requres full mzobjs
local.mz_log_removed_rows <- function( in_mz, out_mz, method){
  mz_before <- nrow(as.data.frame(in_mz))
  mz_after <- nrow(as.data.frame(out_mz))
  mz_removed <- mz_before - mz_after

  print(paste("INFO:", method, ": Before Rows  : ", mz_before, sep=""))
  print(paste("INFO:", method, ": After Rows   : ", mz_after, sep=""))
  print(paste("INFO:", method, ": Dropped Rows : ", mz_removed, sep=""))
}

local.mz_polarity_guesser <- function(input, pos_return = "Positive", neg_return = "Negative"){
  input <- head(input)
  try({input[3,] <- colnames(input)})
  if (sum(grep("[^A-Za-z|^][Nn]eg|[Nn]egative",input)) > 0){
    print("INFO:sumR::polarity_guess: Detected Negative")
    ret <- neg_return
  } else if(sum(grep("[^A-Za-z|^][Pp]os|[Pp]ositive",input)) >0){
    ret <- pos_return
    print("INFO:sumR::polarity_guess: Detected Positive")
  } else {
    stop("ERROR: sumR::polarity_guess: Could not guess positive or negative from colnames")
  }
  ret
}

local.save_plot <- function(plot_name, output_dir = paste("output", local.dir_sep() ,Sys.Date(), "")){
  ggsave(paste(output_dir,local.dir_sep(),plot_name,".png",sep = ""))
}

local.export_thread_env <- function(cores, func_name){
  func_name <- gsub("\\(.*\\)","",func_name)
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, varlist = ls(environment(func_name)),
                                envir = environment(func_name))
  } else {
    cl <- NULL
  }
}
