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
#'      "postbin_method" = max,
#'      "tolerance" = 5e-6,
#'      "timeSquash_method" = mean,
#'      "missingness_threshold" = .1,
#'      "intensity_threshold" = 1000,
#'      "target_list" = c()
#'
#' @param files \cr
#'   String: File or directory of mzmL files
#'   List: of mzML files (i.e. list.files(".",pattern = "qc_.*mzML"))
#' @param parallel
#'   Integer: Can I has multithreading? (Need parallel)
#' @param params
#'   List: override any subset (or all) default parameters
#' @returns MzObj
#' @export
mzml_extract_magic <- function(files = getwd(), polarity = c(0,1), params = NULL ){
  start <- Sys.time()
  files <- extract.file_lister(files)
  #TODO: Fix threaded db connections
  #if (parallel && length(polarity) > 1){
    #cl <- local.export_thread_env(2)
  #} else {
    cl <- NULL
  #}
  #tryCatch({
    df <- pblapply( polarity, cl=cl, function(pol){
      pol_eng <- extract.pol_english(pol)
      params <- extract.fill_defaults(params, pol_eng)
      print(paste0("INFO: mzml_extract_magic : Extraction of [",pol_eng,"], beginning"))
      pbapply::pblapply(files, function(file) mzml_extract_file(file, polarity=pol, cl = NULL, magic=T, params=params))
      print(paste0("INFO: mzml_extract_magic : Extraction of [",pol_eng,"], complete"))
      tables <- DBI::dbGetQuery(params$dbconn, "SELECT table_name FROM duckdb_tables;")
      if (! paste0("PREFINAL_",pol_eng) %in% tables[,]) {
        query <- paste0("CREATE TABLE ",paste0("PREFINAL_",pol_eng), " AS (SELECT * FROM ", tables[1,], paste(" FULL OUTER JOIN", tables[2:length(tables[,]),], collapse = " USING (mz)"), " USING (mz) ORDER BY mz)")
        print(paste0("INFO: mzml_extract_magic : creating initial MZ_OBJ of [",pol_eng,"]"))
        DBI::dbExecute(params$dbconn, query)
      } else {
        print(paste0("INFO: mzml_extract_magic : initial MZ_OBJ of [",pol_eng,"], already exists"))
      }
      print(paste0("INFO: mzml_extract_magic : Final Binning [",pol_eng,"], starting"))
      gc()
      out <- DBI::dbGetQuery(params$dbconn, paste0("select * from PREFINAL_", pol_eng))
      #out <- extract.mz_format(out)
      out <- mzT_binning(out, params$postbin_method, paste0('extract_magic_',pol_eng))
      duckdb::dbWriteTable(params$dbconn ,paste0("FINAL_",pol_eng) , out,append = TRUE)
      print(paste0("INFO: mzml_extract_magic : Final Binning [",pol_eng,"], complete"))
      gc()
      #TODO, move "Select * final_" returns outside of this loop. (So neg doesn't bog up memory)
      return(out)
    })
  #TODO: Fix threaded db connections
  #}, finally={
  # 
  #  if (parallel){
  #    local.kill_threads(cl)
  #  }
  #})
  out <- list()
  for( i in 1:length(polarity)) {
    name <- extract.pol_english(polarity[i])
    out[name] <- df[i]
  } 
  print(paste("INFO: mzml_extract_magic SUCCESS, execution complete in:",round(as.numeric(Sys.time()-start, units="mins"),3),"minutes"))
  out
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
mzml_extract_file <- function(file, polarity = "",  magic = T, cl = NULL , params = NULL) {
  pol_eng <- extract.pol_english(polarity)
  params <- extract.fill_defaults(params,pol_eng)
  name <- strsplit(basename(file),"\\.")[[1]][1]
  name <- paste(pol_eng,name,sep="_")
  if( is.null(params$dbconn) || ! name %in% DBI::dbGetQuery(params$dbconn, "show tables")[[1]] ) {
    if (!is.null(cl)) {
      pbo <- pbapply::pboptions(type = "none")
      on.exit(pbapply::pboptions(pbo), add = TRUE)
    }

    #file <- file.path(file)

    if (!file.exists(file)) {
      warning(sprintf("Cannot find file %s", file))
      return(NULL)
    }
    z <- mzR::openMSfile(file)
    meta <- mzR::header(z)
    range <- abs(meta$scanWindowUpperLimit - meta$scanWindowLowerLimit)
    scans <- range >= params$massWindow[1] & range <= params$massWindow[2]
    p <- mzR::peaks(z)
    if ( length(params$target_list) > 0 ){
      p <- lapply(p,function(pi){
          pi_df <- as.data.frame(pi)
          target_bool <- lapply(params$target_list, function(x){abs(pi_df$mz - x)/x < params$target_list_tolerance});
          target_bool <- Reduce("|",target_bool)
          pi[target_bool,]
      })
    }
    if ( !polarity == "" ) {
      #meta$polarity == 1 == Positive
      p <- p[meta["polarity"] == polarity]
      meta <- meta[ meta["polarity"] == polarity,]
    }
    print(sprintf("INFO: Scans found in %s: %i",name, nrow(meta)))
    rts <-lapply(meta$retentionTime, function(x){c("mz",x)})
    l <- pbapply::pblapply(p, cl = cl, centroid.singleScan)
    for( i in seq_along(l) ){ 
      colnames(l[[i]]) <- rts[[i]] 
    }
    out <- extract.mz_format(l)
    if (magic) {
      print(sprintf("INFO: Binning & Filtering %s", name))
      out <- mzT_filtering( out,
                            prebin_method = params$prebin_method,
                            missingness_threshold = params$missingness_threshold,
                            intensity_threshold = params$intensity_threshold,
                            log_name = name)
      print(sprintf("INFO: TimeSquashing %s",name))
      out <- mzT_squash_time(out, timeSquash_method = params$timeSquash_method)
      colnames(out) <- c('mz',name)
    }
    if (class(params$dbconn) == 'duckdb_connection'){
      duckdb::dbWriteTable(params$dbconn ,name , out,append = TRUE)
    } else {
      return(out)
    }
  } else {
    print(paste0("INFO: [", name, "] found in db, skipping"))
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
  if ( !is.null(prebin_method)) {
    mzT <- mzT_binning(mzT,prebin_method,log_name)
  }
  mzT <- mz_filter_low_intensity(mzT,threshold = intensity_threshold, log_name)
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
#' The mzT_squash_time object will align over the scans in a single file and
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
mzT_squash_time <- function(mzT, timeSquash_method = mean, ignore_zeros = T, cl = NULL){
  # This should be handled by filter low intesity
  if( ignore_zeros ){
    mzT[mzT == 0] <- NA
  }
  squash <- pbapply::pbapply(dplyr::select(mzT, -mz),1, timeSquash_method, na.rm=T, cl=cl)
  squash[is.na(squash)] <- 0
  out <- data.table::as.data.table(cbind(mzT$mz, squash))
  names(out) <- c("mz","intensity")
  out
}

##TODO add notes
#' @export
mzT_binning <- function(df, method = 'max', log_name = "", tolernace = 5e-6){
  db_name <- paste0("Binning_", log_name)
  db_file <- paste0(local.output_dir(),'/tmp-',db_name,'.duckdb')
  bin <- binning(df$mz, tolerance = 5e-6)
  local.mz_log_removed_rows((bin), unique(bin), db_name)
  df$mz <- bin
  rm(bin); gc()
  #TODO see if we can eliminate this check (scanTime as string from extract_file gets weird)
  if (regexpr("extract_magic.*",log_name)){
    cols <- colnames(dplyr::select(df, -mz))
  } else {
    cols <- as.numeric(colnames(dplyr::select(df, -mz)))
  }
  cols <- sprintf('%s("%s")', method, paste(cols ,collapse=paste0('"),',method,'("')))
  #TODO: investigate possible duplications
  query <- paste0("CREATE TABLE final_tmp AS (SELECT mz, ",cols," FROM ", db_name, " GROUP BY mz ORDER by mz)" )
  #Create tmp db to perform binning much quicker via SQL (and with limited memory usage)
  dbconn <- DBI::dbConnect(duckdb::duckdb(db_file))
  duckdb::duckdb_register(dbconn, db_name, df, overwrite = T)
  rm(df); gc()
  DBI::dbExecute(dbconn, query)
  df <- DBI::dbGetQuery(dbconn, "SELECT * FROM final_tmp")
  DBI::dbExecute(dbconn, "DROP TABLE final_tmp")
  DBI::dbDisconnect(dbconn)
  unlink(db_file)
  colnames(df)<- gsub(paste0(method,"\\((.*)\\)"),"\\1",colnames(df))
  df[df < 0] <- 0
  df[is.na(df)] <- 0
  df
}

#convert annoying mzR::list[[mz,i],[mz,i],[mz,i]] into usuable dataframe
extract.mz_format <- function(out){
  out <- dplyr::bind_rows(out)
  out <- out[order(out$mz),]
  out$mz <- as.numeric(out$mz)
  out$mz <- round(out$mz, 5)
  data.table::as.data.table(out)
}

extract.file_lister <- function(files) {
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

extract.pol_english <- function(polarity){
  if(polarity == 1){
    return("pos")
  }else if(polarity == 0){
    return("neg")
  }
}

#TODO move to tools get_defaults
extract.fill_defaults <- function(params = NULL, db_prefix = "all"){
  defaults <- list(
    "dbconn" = duckdb::dbConnect(duckdb::duckdb(paste0(local.output_dir(),'/',db_prefix,'_extract.duckdb'))),
    "prebin_method" = "max",
    "postbin_method" = "max",
    "tolerance" = 5e-6,
    "timeSquash_method" = mean,
    "missingness_threshold" = 1,
    "missingness_threshold" = .1,
    "intensity_threshold" = 1000,
    "target_list" = c(),
    "target_list_tolerance" = 5e-5
  )

  out <- append( params, defaults[is.na(is.na(params)[names(defaults)])] )
  if (! out$prebin_method %in% c('max','min','sum','avg') & out$postbin_method %in% c('max','min','sum','avg')) {
    warning("WARN: fill_defaults: Binning method should be SQL aggregate function (max, min, sum, avg)")
  }

  out
}