# ***  -----------------------------------------------------
#' Create MzObj from directory or file
#'   Current Defaults:
#'      mzBinningMath: max
#'      mzBinTolerance = 5e-6
#'      timeSquashMath: mean
#'      missingness_threshold: 0.02 (2%)
#'      massWindow: c(0,Inf)
#' @param files \cr
#'   String: File or directory of mzmL files
#'   List: of mzML files (i.e. list.files(".",pattern = "qc_.*mzML"))
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#' @param verbose Which Parallization method: parLapply is faster (on Win), but no active output)
#'   Boolean: T=pblapply; F=parLapply
#'
#' @export
mzml_extract_magic <- function(files = getwd(), cores = 4, verbose = F, ... ){
  #[0|1]: 0=Negative, 1=Positive
  magic.polarity_loop <- function(polarity){
    if(polarity){
      pol_eng="pos"
    }else{
      pol_eng="neg"
    }
    print(sprintf("INFO: %s.TIVE",toupper(pol_eng)))
    cl <- local.export_thread_env(cores,environment())
    mzT <- pbapply::pblapply(files, mzml_extract_file, cl=cl)
    print(paste("INFO: mzml_extract_magic :", format(Sys.time(),"%H:%M:%S"), ": Extraction of [",pol_eng,"] complete, now formatting data.",sep=""))
    for( i in seq_along(mzT) ){
      #crappy gsub to change col names to [polarity]_[filename(without extension)]
      colnames(mzT[[i]]) <- c("mz", gsub("^(.*)\\.(.*)",paste(pol_eng,"\\1",sep="_"),files[i]))
    }
    mz <- priv.mz_format(mzT)
    mz <- priv.binning(mz, postbin_method)
  }
  
  
  if (class(files) ==  "character" && any(!grepl('\\.mzML$',files))){
    files <- list.files(path=files, pattern="*.mzML")
  }
  if (length(files) == 0) {
    warning("Cannot find mzML files in given files.")
    return(NULL)
  }

  tryCatch({
    # Prepare each file (by polarity) #TODO, maybe you know...don't
    # Get MzT for each file
    if(cores > 1){
      cl <- local.export_thread_env(2, environment())
      if(verbose){
        out <- pbapply::pblapply(c(0,1), magic.polarity_loop)
      } else {
        out <- parallel::parLapply(cl, c(0,1),magic.polarity_loop)
      }
    } else {
      mz_pos <- magic.polarity_loop(1)
      print("INFO: mzml_extract_magic: Positive Complete")
      mz_neg <- magic.polarity_loop(0)
      out <-list(pos=mz_pos, neg=mz_neg)
    }
    return(out)
  },
  finally={
    local.kill_threads(cl)
  })
}

# ***  -----------------------------------------------------
#' Return mzT from a single file
#'
#' @param file \cr
#'   String : Filename.mzML
#' @param magic \cr
#'   Boolean : To run default filtering, and TimeSquashing
#' @param massWindow \cr
#'   c(Int,Int): Filter the mass window
#' @param polarity \cr
#'   [0|1]: 0=Negative, 1=Positive, NULL=both
#' @param cl \cr
#'   parallel::threadCluster (optional)
#'
#' @export
mzml_extract_file <- function(file, polarity = NULL,  magic = T, massWindow = c(0, Inf),  cl = NULL ) {

  .get_mzt <- function(data, meta, magic = T, cl = cl) {
    #TODO, get this to output 'file' correctly
    #TODO, unnest this. I keep getting this error:
    #  cannot coerce type 'closure' to vector of type 'character'
    print(sprintf("INFO: %s: Scans found in %s: %i",format(Sys.time(),"%H:%M:%S"), file, nrow(meta)))
    rts <-lapply(meta$retentionTime, function(x){c("mz",x)})
    l <- pbapply::pblapply(data, cl = cl, centroid.singleScan)
    for( i in seq_along(l) ){ colnames(l[[i]]) <- rts[[i]] }
    out <- priv.mz_format(l)
    if (magic) {
      out <- sumR::mzT_filtering(out)
      out <- mzT_squashTime(out, cl=cl)
    }else{
      out
    }
  }

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
  df <- mzR::header(z)
  range <- abs(df$scanWindowUpperLimit - df$scanWindowLowerLimit)
  scans <- range >= massWindow[1] & range <= massWindow[2]
  p <- mzR::peaks(z)
  if ( !is.null(polarity) ) {
    #df$polarity == 1 == Positive

    out <- .get_mzt(p[df["polarity"] == polarity],
                    df[ df["polarity"] == polarity,],
                    magic = magic, cl = cl)
  } else {
    out <- .get_mzt(p, df, cl = cl)
  }
}

# ***  -----------------------------------------------------
#' MzTime filtering (Uses mz_filtering) & binning for an MzT obj.
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
#' @export
# Note missingness_threshold is very low
mzT_filtering <- function(mzT, prebin_method = max, missingness_threshold = 1, intensity_threshold = 1000 ){
  mzT <- priv.binning(mzT,prebin_method)
  mzT <- mz_filter_lowIntensity(mzT,threshold = intensity_threshold)
  mzT <- mz_filter_missingness(mzT,threshold = missingness_threshold)
}

# ***  -----------------------------------------------------
#' Squash MzT obj (many scans) into a single mz|intensity dataframe
#'
#' @param mzT \cr
#'   MzT : MzT object
#' @param method \cr
#'   [Math] : i.e. (mean, max, median)
#' @param ignore_zeros \cr
#'   Boolean: Should we set 0 <- NA (to not affect the math)
#'
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

priv.binning <- function(df, method = max, tolerance = 5e-6){
  bin <- binning(df$mz, tolerance = tolerance)
  local.mz_log_removed_rows((bin), unique(bin), "Binning")
  df <- suppressWarnings(aggregate(dplyr::select(df,-mz),list(bin), method, na.rm = T, na.action=NULL))
  names(df)[1] <- "mz"
  df[df < 0] <- 0
  df[is.na(df)] <- 0
  df
}

#convert annoying mzR::list[[mz,i],[mz,i],[mz,i]] into usuable dataframe
priv.mz_format <- function(out){
  out <- dplyr::bind_rows(out)
  out <- out[order(out$mz),]
  out$mz <- round(out$mz, 5)
  #rownames(out) <- NULL
  out
}
