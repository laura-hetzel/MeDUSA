# ***  -----------------------------------------------------
#' MZML Extract Magic: into MZ_OBJ
#'
#' MZML->MZ_OBJ for all your MZ_OBJ needs
#'  - Requires: Bioc::mzR, dplyr, pbapply
#'
#' @param input_mzml \cr
#'   String : Filename.mzml
#' @param threads \cr
#'   Numeric : How many threads?
#' @returns list[mz_obj-pos, mz_obj-neg]
#'
#' @export
#
mzml_magic <- function(input_dir = getwd(), threads=2){
  .polarity_loop <- function(df,polarity){
    df <- pbapply::pblapply(df, cl=threads, .bin)
    out <- Reduce(function(x,y) merge(x,y, by='mz', all=T),df)
    print(paste("INFO: MZ-OBJ created of dimension: ", dim(out)[1],", ",dim(out)[2],sep=""))
    #TODO, decide if we want to manage neg/pos this way
    names(out) <- c("mz", paste(names(df),polarity,sep="_"))
    out <- .rebin(out)
    out[ is.na(out) ] <- 0
    out
  }
  mzs <- .old_magic(input_dir)
  pos_out <- .polarity_loop(mzs$pos,"pos")
  neg_out <- .polarity_loop(mzs$neg,"neg")
  list(pos = pos_out, neg = neg_out)
}

#Bin each sample individually ###TODO find source
###TODO: Choose bin intesity by (average, max, and weighted-intensity)
.bin <- function(x, math = max){
  x<-dplyr::select(x,mz,i)
  bin <- binning(x$mz)
  x<-aggregate(x$i,list(bin),math,na.rm=T)
  names(x)<-c("mz","intensity")
  x
}

#Re-Bin after the Mz-Obj is created
.rebin <- function(x, math = max){
  bin <- binning(x$mz)
  local.mz_log_removed_rows((bin), unique(bin), "Re-Binning")
  x <- aggregate(dplyr::select(x,-mz),list(bin),math)
  names(x)[1] <- "mz"
  x
}

#Extract, centroid, TimeSquash
.old_magic <- function(input_dir = getwd()){
  files <- list.files(path=input_dir, pattern="*.mzML")
  files <- sapply(files, function(x){paste(input_dir,x, sep="/")})

  #uses centroiding.R; requires: pbapply, tools, mzR, parallel?
  mz_pos <- extractPeaks(files,polarity = "+", combineSpectra=T)
  mz_neg <- extractPeaks(files,polarity = "-", combineSpectra=T)
  list(pos = mz_pos, neg = mz_neg)
}

##TODO add timequash math combineSpectra=F
