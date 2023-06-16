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
#' @param cores
#'   Integer: Can I has multithreading? (Need parallel)
#' @export
mzml_extract_magic <- function(files = getwd(), cores = 1,  ...){
  if (!grepl('\\.mzML$',files)){
    files <- list.files(path=files, pattern="*.mzML")
  }
  if (length(files) == 0) {
    warning("Cannot find mzML files in given files.")
    return(NULL)
  }
  tryCatch({
    # Make cluster if cores > 1
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
    } else {
      cl <- NULL
    }

    # Prepare each file (by polarity) #TODO, maybe you know...don't
    # Get MzT for each file
    mz_pos <- .magic_polarity_loop(files,1,cl)
    mz_neg <- .magic_polarity_loop(files,0,cl)
    list(pos=mz_pos, neg=mz_neg)
  },
  finally={
    #This doesn't seem to stop the cluster :/
    if (cores > 1) {
      parallel::stopCluster(cl)
    }
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
    print(sprintf("INFO: Scans found in %s: %i",file, nrow(meta)))
    rts <-lapply(meta$retentionTime, function(x){c("mz",x)})
    l <- pbapply::pblapply(data, cl = cl, centroid.singleScan)
    for( i in seq_along(l) ){ colnames(l[[i]]) <- rts[[i]] }
    out <- .mz_format(l)
    if (magic) {
      out <- mzT_filtering(out)
      out <- mzT_squashTime(out)
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
    data <-
    meta <-

    out <- .get_mzt(p[df["polarity"] == polarity],
                    df[ df["polarity"] == polarity,],
                    magic = magic, cl = cl)
  } else {
    out <- .get_mzt(p, df, cl = cl)
  }
}

# ***  -----------------------------------------------------
#' Missingness filtering (Uses mz_filtering) & binning for an MzT obj.
#'
#' @param mzT \cr
#'   MzT : MzT object
#' @param method \cr
#'   [Math] : i.e. (mean, max, median)
#' @param missingness_threshold \cr
#'   Float: Percent as a decimal. i.e. requiring non-zero values in 0.02 (2%) of columns per mz
#'
#' @export
mzT_filtering <- function(mzT, method = max, missingness_threshold = 0.02){
  #MzT object (pre)binning
  mzT <- .binning(mzT,method)
  mzT <- mz_filter_missingness(mzT,threshold = missingness_threshold)
}

# ***  -----------------------------------------------------
#' Squash MzT obj (many scans) into a single mz|intensity dataframe
#'
#' @param mzT \cr
#'   MzT : MzT object
#' @param method \cr
#'   [Math] : i.e. (mean, max, median)
#' @param ignore_nulls \cr
#'   Boolean: Consider 0 == null, thena ignore them?
#'
#' @export
mzT_squashTime <- function(mzT, method = mean, ignore_nulls = T){
  if(ignore_nulls){
    mzT[mzT==0] <- NA
  }
  squash <- apply(dplyr::select(mzT, -mz),1, method, na.rm=ignore_nulls)
  squash[is.na(squash)] <- 0
  out <- as.data.frame(cbind(mzT$mz, squash))
  names(out) <- c("mz","intensity")
  out
}

.magic_polarity_loop <- function(files,pol,cl){
  if(pol){
    pol_eng="pos"
  }else{
    pol_eng="neg"
  }
  print(sprintf("INFO: %s.TIVE",toupper(pol_eng)))
  mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=pol, cl = cl))

  for( i in seq_along(mzT) ){
    #crappy gsub to change col names to [polarity]_[filename(without extension)]
    colnames(mzT[[i]]) <- c("mz", gsub("^(.*)\\.(.*)",paste(pol_eng,"\\1",sep="_"),files[i]))
  }
  mz <- .mz_format(mzT)
  mz <- .binning(mz)
}

.binning <- function(df, method = max, tolerance = 5e-6){
  bin <- binning(df$mz, tolerance = tolerance)
  local.mz_log_removed_rows((bin), unique(bin), "Binning")
  df <- aggregate(dplyr::select(df,-mz),list(bin), method)
  names(df)[1] <- "mz"
  df[is.na(df)] <- 0
  df
}

#convert annoying mzR::list[[mz,i],[mz,i],[mz,i]] into usuable dataframe
.mz_format <- function(out){
  out <- dplyr::bind_rows(out)
  out <- out[order(out$mz),]
  out$mz <- round(out$mz, 5)
  rownames(out) <- NULL
  out
}
