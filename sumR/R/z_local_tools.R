### z_* written by ehetzel zz_* written by others

local.dir_sep <- function() {
  if ( .Platform$OS.type == "unix" ){
    dir_sep = "/" #Unix
  } else {
    dir_sep = "\\" #Windows, note escape char
  }
  dir_sep
}

local.output_dir <- function(){
  paste0("output", local.dir_sep() ,Sys.Date())
}

#Show zeros requres full mzobjs
local.mz_log_removed_rows <- function( in_mz, out_mz, method){
  mz_before <- nrow(as.data.frame(in_mz))
  mz_after <- nrow(as.data.frame(out_mz))
  mz_removed <- mz_before - mz_after

  print(paste0("INFO:", method, ": Before Rows  : ", mz_before))
  print(paste0("INFO:", method, ": After Rows   : ", mz_after))
  print(paste0("INFO:", method, ": Dropped Rows : ", mz_removed))
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

local.meta_polarity_fixer <- function(input_mz, meta){
  prepend <- local.mz_polarity_guesser(input_mz, pos_return = "pos", neg_return = "neg")
  meta$sample_name <- paste(prepend, meta$sample_name, sep="_")
  meta
}

local.save_plot <- function(plot_name, output_dir = local.output_dir(), dim=c(8,8)){
  ggplot2::ggsave(filename = paste0(output_dir,local.dir_sep(),plot_name,".png"),
                  height = dim[1],
                  width = dim[2]
                  )
}

local.export_thread_env <- function(cores, env = environment()){
  if (cores > 1) {
    cl <- parallel::makeCluster(cores, outfile="")
    parallel::clusterExport(cl, varlist = ls(env),
                                envir = env)
  } else {
    cl <- NULL
  }
  return(cl)
}

local.kill_threads <- function(cl = NULL){
  tryCatch(parallel::stopCluster(cl),
           error = function(e){print(paste("WARN: Thread kill error:",e))},
           finally = {print("INFO: Killing threads. Often throws odd errors.")
                     suppressWarnings(try(closeAllConnections(),silent=TRUE))})
}

local.ensure_mz <- function(input_a, input_b, source){
  if (nrow(input_a) != nrow(input_b)){
    stop(paste0("ERROR: ", source, ": Dataframes have a different number of mz"))
  }
  if ( sum(c(colnames(input_a), colnames(input_b)) == "mz")) {
    if ( sum(colnames(input_b) == "mz") ) {
      mz <- input_b$mz
      input_b <- dplyr::select(input_b, -mz)
    }
    if (sum(colnames(input_a) == "mz") ) {
      mz <- input_a$mz
      input_a <- dplyr::select(input_a, -mz)
    }
  } else {
    stop(paste0("ERROR: ", source, ": neither input dataframes have column mz"))
  }
  list(mz   = as.data.frame(mz), 
       df_a = as.data.frame(input_a), 
       df_b = as.data.frame(input_b)) 
}


