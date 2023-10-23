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
#'
#' @export
mzml_extract_magic <- function(files = getwd(), cores = 6,  params = NULL ){
  start <- Sys.time()
  params <- magic.fill_defaults(params)

  if (class(files) ==  "character" & any(!grepl('\\.mzML$',files))){
    files <- list.files(path=files, pattern="*.mzML")
  }
  if (length(files) == 0) {
    warning("Cannot find mzML files in given files.")
    return(NULL)
  }

  tryCatch({
    if (cores > 1){
      cl <- local.export_thread_env(2)
      out <- parallel::parLapply(cl, c(0,1), function(x) magic.polarity_loop(polarity = x, files, cores, params))
      names(out) <- c('neg','pos')
    } else {
      mz_pos <- magic.polarity_loop(files,1,params)
      print("INFO: mzml_extract_magic: Positive Complete")
      mz_neg <- magic.polarity_loop(files,0,params)
      out <- list(pos=mz_pos, neg=mz_neg)
    }
    print(paste("INFO: mzml_extract_magic, execution complete in:",round(as.numeric(Sys.time()-a, units="mins"),3),"minutes"))
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
#' @param polarity \cr
#'   [0|1]: 0=Negative, 1=Positive, NULL=both
#' @param cl \cr
#'   parallel::threadCluster (optional)
#'
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
mzT_filtering <- function(mzT, prebin_method = max, missingness_threshold = 1, intensity_threshold = 1000, log_name = "" ){
  mzT <- magic.binning(mzT,prebin_method,log_name)
  mzT <- mz_filter_lowIntensity(mzT,threshold = intensity_threshold, log_name)
  mzT <- mz_filter_missingness(mzT,threshold = missingness_threshold,log_name)
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
mzT_squashTime <- function(mzT, timeSquash_method = mean, ignore_zeros = T){
  # This should be handled by filter low intesity
  if( ignore_zeros ){
    mzT[mzT == 0] <- NA
  }
  squash <- apply(dplyr::select(mzT, -mz),1, timeSquash_method, na.rm=T)
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
  if(cores > 3 ){
    cl <- local.export_thread_env(round(cores/2) , environment(magic.polarity_loop))
    mzT <- parallel::parLapply(cl, files, function(file) mzml_extract_file(file, polarity, T,  NULL , params))
  } else {
    mzT <- pbapply::pblapply(files,  function(x) mzml_extract_file(x, polarity, cl))
  }
  print(paste("INFO: mzml_extract_magic : Extraction of [",pol_eng,"] complete, now formatting data.",sep=""))
  for( i in seq_along(mzT) ){
    #crappy gsub to change col names to [polarity]_[filename(without extension)]
    colnames(mzT[[i]]) <- c("mz", gsub("^(.*)\\.(.*)",paste(pol_eng,"\\1",sep="_"),files[i]))
  }
  mz <- magic.mz_format(mzT)
  mz <- magic.binning(mz, params$postbin_method, "final")
}

magic.binning <- function(df, method = max, log_name = "", tolerance = 5e-6){
  bin <- binning(df$mz, tolerance = tolerance)
  local.mz_log_removed_rows((bin), unique(bin), paste("Binning",log_name))
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
  #rownames(out) <- NULL
  out
}

magic.fill_defaults <- function(params){
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
