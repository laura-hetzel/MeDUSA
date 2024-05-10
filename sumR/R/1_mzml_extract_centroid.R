# ***  -----------------------------------------------------
#' Create MzObj from directory or file (with filtering & binning)
#'
#' @description
#' The mzml file contains data that is not aligned. Two alignments are
#' necessary for a proper comparison of data.
#' The first alignment compares data in one file, over multiple scans. This
#' alignment removes the time factor and creates one intensity per m/z value for
#' each file. The default is to create a mean of the intensity per m/z. A
#' missingness threshold is set so that an m/z must occur in a set percentage of
#' scans in order to not be filtered out as an artifact. The tolerance for an
#' m/z to be considered the same feature across multiple scans is referred to as
#' the mzBinTolerance. If m/z values within the bin differ, the mzBinningMath
#' will determine which method is used to determine the assigned m/z value.
#' The second alignment compares the data in all of the files being processed.
#' The same steps and math are used to determine the m/z values for all files.
#' The intensity does not undergo further math and is simply assigned to the
#' corresponding m/z value.
#' The mzml_extract_magic function performs the above mentioned alignment in one
#' function.
#'
#'   Current Defaults:
#'      mzBinningMath: max
#'      mzBinTolerance = 5e-6
#'      timeSquashMath: mean
#'      missingness_threshold: 0.02 (2%)
#'      massWindow: c(0,Inf)
#'
#' @param files \cr
#'   String: File or directory of mzmL files
#'   List: of mzML files (i.e. list.files(".",pattern = "qc_.*mzML"))
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#' @param verbose Which Parallization method: parLapply is faster (on Win), but no active output)
#'   Boolean: T=pblapply; F=parLapply
#' @returns MzObj
#' @export
mzml_extract_magic <- function(files = getwd(), cores = 6,  params = NULL ){
  start <- Sys.time()
  params <- magic.fill_defaults(params)
  files <- magic.file_lister(files)

  tryCatch({
    if (cores > 1){
      cl <- local.export_thread_env(2)
      out <- parallel::parLapply(cl, c(0,1), function(x) magic.polarity_loop(polarity = x, files, cores, params))
      names(out) <- c('neg','pos')
      local.kill_threads(cl)
      gc()
      cl <- local.export_thread_env(2)
      out <- parallel::parLapply(cl, c("neg","pos"), function(x) magic.polarity_final(x, out, files, params))
      names(out) <- c('neg','pos')
    } else {
      mz_pos <- magic.polarity_loop(files,1,params)
      print("INFO: mzml_extract_magic: Positive Complete")
      mz_neg <- magic.polarity_loop(files,0,params)
      out <- list(pos=mz_pos, neg=mz_neg)
    }
    print(paste("INFO: mzml_extract_magic SUCCESS, execution complete in:",round(as.numeric(Sys.time()-start, units="mins"),3),"minutes"))
    return(out)
  },
  finally={
    local.kill_threads()
  })
}

# ***  -----------------------------------------------------
#' Return MzT (mzTime Object) from a single file
#'
#' @description
#' The mzml file contains data that is not aligned. Two alignments are
#' necessary for a proper comparison of data.
#' The first alignment compares data in one file, over multiple scans. This
#' alignment removes the time factor and creates one intensity per m/z value for
#' each file. The default is to create a mean of the intensity per m/z. A
#' missingness threshold is set so that an m/z must occur in a set percentage of
#' scans in order to not be filtered out as an artifact. The tolerance for an
#' m/z to be considered the same feature across multiple scans is referred to as
#' the mzBinTolerance. If m/z values within the bin differ, the mzBinningMath
#' will determine which method is used to determine the assigned m/z value.
#' The second alignment compares the data in all of the files being processed.
#' The same steps and math are used to determine the m/z values for all files.
#' The intensity does not undergo further math and is simply assigned to the
#' corresponding m/z value.
#' The mzml_extract_file function will only extract the data from the mzml files
#' and create a dataframe that is compatible with the other functions in this
#' package. No alignment or other processing is performed on the data.
#'
#' @param file \cr
#'   String : Filename.mzML
#' @param magic \cr
#'   Boolean : To run default filtering, and TimeSquashing
#' @param polarity \cr
#'   [0|1]: 0=Negative, 1=Positive, NULL=both
#' @param cl \cr
#'   parallel::threadCluster (optional)
#' @returns MzT object
#' @export
mzml_extract_file <- function(file, polarity = "",  magic = T,  cl = NULL, params = NULL) {
  params <- magic.fill_defaults(params)

  if (!is.null(cl)) {
    pbo <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(pbo), add = TRUE)
  }

  file <- file.path(file)

  if (!file.exists(file)) {
    warning(sprintf("Cannot find file %s", file))
    return(NULL)
  }
  z <- mzR::openMSfile(file)
  meta <- mzR::header(z)
  range <- abs(meta$scanWindowUpperLimit - meta$scanWindowLowerLimit)
  scans <- range >= params$massWindow[1] & range <= params$massWindow[2]
  p <- mzR::peaks(z)
  if ( !polarity == "" ) {
    #meta$polarity == 1 == Positive

    p <- p[meta["polarity"] == polarity]
    meta <- meta[ meta["polarity"] == polarity,]

    if(polarity){
      polarity="pos"
    }else{
      polarity="neg"
    }
  }

  print(sprintf("INFO: Scans found in %s %s: %i",file, polarity, nrow(meta)))
  rts <-lapply(meta$retentionTime, function(x){c("mz",x)})
  l <- pbapply::pblapply(p, cl = cl, centroid.singleScan)
  for( i in seq_along(l) ){ colnames(l[[i]]) <- rts[[i]] }
  out <- magic.mz_format(l)
  if (magic) {
    out <- mzT_filtering( out,
                          prebin_method = params$prebin_method,
                          missingness_threshold = params$missingness_threshold,
                          intensity_threshold = params$intensity_threshold,
                          log_name = paste(file,polarity))
    out <- mzT_squashTime(out, timeSquash_method = params$timeSquash_method)
  }else{
    out
  }
}

# ***  -----------------------------------------------------
#' MzTime filtering (Uses mz_filtering) & binning for an MzT obj.
#'
#' @description
#' The mzml file contains data that is not aligned. Two alignments are
#' necessary for a proper comparison of data.
#' The first alignment compares data in one file, over multiple scans. This
#' alignment removes the time factor and creates one intensity per m/z value for
#' each file. The default is to create a mean of the intensity per m/z. A
#' missingness threshold is set so that an m/z must occur in a set percentage of
#' scans in order to not be filtered out as an artifact. The tolerance for an
#' m/z to be considered the same feature across multiple scans is referred to as
#' the mzBinTolerance. If m/z values within the bin differ, the mzBinningMath
#' will determine which method is used to determine the assigned m/z value.
#' The second alignment compares the data in all of the files being processed.
#' The same steps and math are used to determine the m/z values for all files.
#' The intensity does not undergo further math and is simply assigned to the
#' corresponding m/z value.
#' The mzT_filtering function will filter the data that is produced removing
#' m/z values that do not fulfill the missing requirement or the minimum
#' intensity threshold.
#'
#'
#' @param mzT \cr
#'   MzT : MzT object
#' @param method \cr
#'   [Math] : i.e. (mean, max, median)
#' @param missingness_threshold \cr
#'   Numeric   : If Decimal : Threshold % to blank data under. (as a decimal 10% = 0.1)
#'             : If Whole : Number of scans required to be nonzero ( Value required in at least 2 scans )
#' @param intensity_threshold \cr
#'   Numeric   : Lowest allowed intensity
#' @returns MzT object
#' @export
# Note missingness_threshold is very low
mzT_filtering <- function(mzT, prebin_method = max, missingness_threshold = 1, intensity_threshold = 1000, log_name = "" ){
  mzT <- magic.binning(mzT,prebin_method,log_name)
  mzT <- mz_filter_lowIntensity(mzT,threshold = intensity_threshold, log_name)
  mzT <- mz_filter_missingness(mzT,threshold = missingness_threshold,log_name)
}

# ***  -----------------------------------------------------
#' "Squash" MzT scans into a single MzObj Intensity
#'
#' @description
#' The mzml file contains data that is not aligned. Two alignments are
#' necessary for a proper comparison of data.
#' The first alignment compares data in one file, over multiple scans. This
#' alignment removes the time factor and creates one intensity per m/z value for
#' each file. The default is to create a mean of the intensity per m/z. A
#' missingness threshold is set so that an m/z must occur in a set percentage of
#' scans in order to not be filtered out as an artifact. The tolerance for an
#' m/z to be considered the same feature across multiple scans is referred to as
#' the mzBinTolerance. If m/z values within the bin differ, the mzBinningMath
#' will determine which method is used to determine the assigned m/z value.
#' The second alignment compares the data in all of the files being processed.
#' The same steps and math are used to determine the m/z values for all files.
#' The intensity does not undergo further math and is simply assigned to the
#' corresponding m/z value.
#' The mzT_squashTime object will align over the scans in a single file and
#' produce one intensity per m/z value per sample.
#' Squash MzT obj (many scans) into a single mz|intensity dataframe
#'
#' @param mzT \cr
#'   MzT : MzT object
#' @param method \cr
#'   [Math] : i.e. (mean, max, median)
#' @param ignore_zeros \cr
#'   Boolean: Should we set 0 <- NA (to not affect the math)
#' @returns MzObj of one sample column
#' @export
mzT_squashTime <- function(mzT, timeSquash_method = mean, ignore_zeros = T, cl = NULL){
  # This should be handled by filter low intesity
  if( ignore_zeros ){
    mzT[mzT == 0] <- NA
  }
  squash <- pbapply::pbapply(dplyr::select(mzT, -mz),1, timeSquash_method, na.rm=T, cl=cl)
  squash[is.na(squash)] <- 0
  out <- as.data.frame(cbind(mzT$mz, squash))
  names(out) <- c("mz","intensity")
  out
}

#[0|1]: 0=Negative, 1=Positive
magic.polarity_loop <- function(files, polarity, cores, params){
  if(polarity){
    pol_eng="pos"
  }else{
    pol_eng="neg"
  }
  print(sprintf("INFO: %s.TIVE",toupper(pol_eng)))
  if(cores > 3 && polarity ){
    #Give Pos an extra core
    cl <- trunc(cores/2)+1
    if(!polarity){
      cl <- cores- cl
    }
    print(paste("cores",polarity,cl))
    cl <- local.export_thread_env(cl , environment(magic.polarity_loop))
    mzT <- parallel::parLapplyLB(cl, files, function(file) mzml_extract_file(file, polarity, T,  NULL , params))
  } else {
    mzT <- pbapply::pblapply(files,  function(x) mzml_extract_file(x, polarity, T, NULL, params))
  }
  print(paste0("INFO: mzml_extract_magic : Extraction of [",pol_eng,"] complete, now formatting data."))
  #Useful for debugging
  #save(mzT, file = paste(Sys.Date(),"-mzT-preBin",pol_eng,".Rdata", sep=""))
  mzT
}

magic.polarity_final <- function(pol,mzT,files,params){
  mzT <- mzT[[pol]]
  for( i in seq_along(mzT) ){
    #crappy gsub to change col names to [polarity]_[filename(without extension)]
    #  Needs to be a for loop, because all the "intensity" columns have to be; lest they be overwritten by mz_format
    colnames(mzT[[i]]) <- c("mz", gsub("(^.*/|^)(.*).mzML",paste(pol,"\\2",sep="_"),files[i]))
  }
  mz <- magic.mz_format(mzT)
  mz <- magic.binning(mz, params$postbin_method)
  #colnames(mz) <- c("mz", gsub("^(.*)\\.(.*)",paste(pol,"\\1",sep="_"),files))
}

magic.binning <- function(df, method = max, log_name = "", tolerance = 5e-6){
  bin <- binning(df$mz, tolerance = tolerance)
  #local.mz_log_removed_rows((bin), unique(bin), paste("Binning",log_name))
  df <- suppressWarnings(aggregate(dplyr::select(df,-mz),list(bin), method, na.rm = T, na.action=NULL))
  names(df)[1] <- "mz"
  df[df < 0] <- 0
  df[is.na(df)] <- 0
  df
}

#convert annoying mzR::list[[mz,i],[mz,i],[mz,i]] into usuable dataframe
magic.mz_format <- function(out){
  out <- dplyr::bind_rows(out)
  out <- out[order(out$mz),]
  out$mz <- round(out$mz, 5)
  as.data.frame(out)
}

magic.file_lister <- function(files) {
  #if string, and not an mzml file, treat as directory
  if (class(files) ==  "character" & any(!grepl('\\.mzML$',files))){
    files <- list.files(path=files, pattern="*.mzML", full.names=T)
  }
  if (length(files) == 0) {
    warning("ERROR: mzml_extract_magic: Cannot find mzML files in given files.")
    return(NULL)
  }
  files
}

magic.fill_defaults <- function(params = NULL){
  defaults <- list(
    "prebin_method" = max,
    "postbin_method" = max,
    "tolerance" = 5e-6,
    "timeSquash_method" = mean,
    "missingness_threshold" = .1,
    "intensity_threshold" = 1000,
    "massWindow" = c(0, Inf)
  )
  append( params, defaults[is.na(is.na(params)[names(defaults)])] )
}
