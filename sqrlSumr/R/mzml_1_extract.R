# ***  -----------------------------------------------------
#' MZML Extract into mz-T(time) OBJ
#' mzT OBJ= dataframe with an mz column, and a column for each scan time's intensity.
#'
#' Extracts MZML file, and splits them into pos/neg scan lists
#'  - Requires: Bioc::mzR
#'
#' @param input_mzml \cr
#'   String : Filename
#' @returns [mzT-pos, mzT-neg]
#'
#' @export
mzml_extract <- function(input_mzml) {
  mzML_xml <- mzR::openMSfile(input_mzml,  backend = "pwiz")
  mzML_header <- mzR::header(mzML_xml)

  print(paste("INFO: converting Pos-file:", input_mzml))
  mzML_posScans <- mzML_header[mzML_header$polarity == 1,]$seqNum
  mzML_posList <- mzR::spectra(mzML_xml,mzML_posScans)
  #what is retention vs injection time?
  names(mzML_posList) <- mzML_header[mzML_header$polarity == 1,]$retentionTime
  mzT_pos <- .mzT_convert(mzML_posList)

  #Annoying splitting neg/pos for memory management
  print(paste("INFO: converting Neg-file:", input_mzml))
  mzML_negScans <- mzML_header[mzML_header$polarity == 0,]$seqNum
  mzML_negList <- mzR::spectra(mzML_xml,mzML_negScans)
  names(mzML_negList) <- mzML_header[mzML_header$polarity == 0,]$retentionTime
  mzT_neg <- .mzT_convert(mzML_negList)

  list(pos = mzT_pos,neg = mzT_neg)
}

# ***  -----------------------------------------------------
#' MzT OBJ Squash mz (simple binning)
#'
#' Squash the mz element of MzT OBJ.
#'  - Requires: Dplyr?
#'
#' @param mzT \cr
#'   mzT Object
#' @param method
#'   String: Mean? Median? etc?
#' @param bin_size \cr
#'   Float: Allowed size of bins
#' @param ignore_nulls \cr
#'   Boolean: Ignore null in the calculation?
#' @returns [mzT-pos, mzT-neg]
#'
#' @export
mzTime_squashMz <- function(mzT, method = mean, bin_size = 0.005, ignore_nulls = TRUE){
  precision <- nchar(sub('.*\\.', '', bin_size))
  if ( max(mzT$mz)-min(mzT$mz) < bin_size ){
    stop("ERROR: bin_size must be greater than the difference of mz")
  }
  if(ignore_nulls){
    mzT[mzT==0] <- NA
  }
  squash <- aggregate(mzT, by=list(cut(mzT$mz , seq(min(mzT$mz), max(mzT$mz), by = bin_size))), method, na.rm=ignore_nulls)
  squash$mz <- round(squash$mz,precision)
  squash[is.na(squash)] <- 0
  names(squash) <- sub("\\..*","",names(squash))
  dplyr::select(squash, -Group)
}

# ***  -----------------------------------------------------
#' MzT OBJ Squash time
#'
#' Squash the time element of MzT OBJ.
#'  - Requires: Dplyr
#'
#' @param mzT \cr
#'   mzT Object
#' @param method
#'   String: Mean? Median? etc?
#' @param ignore_nulls \cr
#'   Boolean: Ignore null in the calculation?
#' @returns [mzT-pos, mzT-neg]
#'
#' @export
mzTime_squashTime <- function(mzT, method = mean, ignore_nulls = TRUE) {
  if(ignore_nulls){
    mzT[mzT==0] <- NA
  }
  squash <- apply(dplyr::select(mzT, -mz),1, method, na.rm=ignore_nulls)
  squash[is.na(squash)] <- 0
  out <- as.data.frame(cbind(mzT$mz, squash))
  names(out) <- c("mz","intensity")
  out
}

# ***  -----------------------------------------------------
#' MZML Extract Magic: into MZ_OBJ
#'
#' MZML->MZ_OBJ for all your MZ_OBJ needs
#'  - Requires: Bioc::mzR, dplyr
#'
#' @param input_mzml \cr
#'   String : Filename.mzml
#' @returns list[mz_obj-pos, mz_obj-neg]
#'
#' @export
#
mzml_magic <- function(input_dir = getwd()){
  .add_mz_col <- function(file, list, polarity = 'pos'){
    local_out <- mzTime_squashTime(list[[polarity]])
    names(local_out) <- c("mz",file)
    output <- readRDS(paste("tmp_mz_",polarity,".RDS", sep=""))
    output <- dplyr::full_join(output, local_out, by="mz")
    saveRDS(mz_pos,file="tmp_mz_pos.RDS")
  }

  mz_pos <- data.frame(mz = double())
  mz_neg <- data.frame(mz = double())
  saveRDS(mz_pos,file="tmp_mz_pos.RDS")
  saveRDS(mz_neg,file="tmp_mz_neg.RDS")
  files <- list.files(path=input_dir, pattern="*.mzML")

  for(file in files){
    mz_list <- mzml_extract(paste(input_dir,file, sep="/"))
    .add_mz_col(file, mz_list, "pos")
    gc()
    .add_mz_col(file, mz_list, "neg")
    gc()
  }

  mz_pos <- readRDS("tmp_mz_pos.RDS")
  mz_neg <- readRDS("tmp_mz_neg.RDS")
  pos_out <- mzTime_squashMz(pos_out)
  neg_out <- mzTime_squashMz(neg_out)
  #file.remove('tmp_mz_pos.RDS','tmp_mz_neg.RDS')
  list(pos = pos_out, neg = neg_out)
}


.mzT_convert <- function(spectra_list) {
  output <- data.frame(mz = double())
  for( columnName in names(spectra_list) ) {
    out <- as.data.frame(spectra_list[[as.character(columnName)]])
    names(out) <- c('mz',as.character(columnName))
    out$mz <- round(out$mz,5)
    output <- dplyr::full_join(output, out, by="mz")
    print(paste("INFO: converting scanTime:", columnName))
    #print(paste("INFO: Current MZ_T Size", dim(output)))
    #TODO I hate this, find a better way to manage memory
    gc()
  }
  output[is.na(output)] <- 0
  output
}
##I'd rather use sapply, but i dont think join like it
#.df_from_list <- function(columnName)
#These joins are expensive; perhaps prebin?
#output <- sapply(names(spectra_list)), .df_from_list)
