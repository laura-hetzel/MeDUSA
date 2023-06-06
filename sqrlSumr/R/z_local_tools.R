### z_* written by ehetzel zz_* written by others

dir_sep = "/" #Unix
#dir_sep = "\" #Windows

local.tmp_dir_exists <- function(name) {
 file.exists(paste(Sys.getenv("tmp_dir"),name, sep = Sys.getenv("dir_sep")))
}

local.tmp_dir_write <- function(filename,sheet,value){
  Tmp_dir <- paste(Sys.getenv("tmp_dir"),name, sep = Sys.getenv("dir_sep"))
  write_excel(value, file = filename, sheet = sheet )
}

local.tmp_dir_read <- function(name, sheet){
  Tmp_dir <- paste(Sys.getenv("tmp_dir"),name, sep = Sys.getenv("dir_sep"))
  data.frame(read_excel(name, sheet = sheet, col_names = T))
}

local.mz_log_removed_rows <- function(in_mz, out_mz, method){
  mz_before <- nrow(as.data.frame(in_mz))
  mz_after <- nrow(as.data.frame(out_mz))
  mz_removed <-  mz_before - mz_after
  print(paste("INFO:", method ,": Before Rows  : ", mz_before, sep=""))
  print(paste("INFO:", method ,": After Rows   :  ", mz_after, sep=""))
  print(paste("INFO:", method ,": Dropped Rows : ", mz_removed, sep=""))
}

local.mz_polarity_guesser <- function(input, pos_return = "Positive", neg_return = "Negative"){
  input <- head(input)
  try({input[3,] <- colnames(input)})
  if (sum(grep("[^A-Za-z|^][Nn]eg|[Nn]egative",input)) > 0){
    print("INFO:sqrlSumr::polarity_guess: Detected Negative")
    ret <- neg_return
  } else if(sum(grep("[^A-Za-z|^][Pp]os|[Pp]ositive",input)) >0){
    ret <- pos_return
    print("INFO:sqrlSumr::polarity_guess: Detected Positive")
  } else {
    stop("ERROR: sqrlSumr::polarity_guess: Could not guess positive or negative from colnames")
  }
  ret
}

local.save_plot <- function(plot_name, output_dir = paste("output",dir_sep,Sys.Date(),"")){
  ggsave(paste(output_dir,dir_sep,plot_name,".png",sep = ""))
}
