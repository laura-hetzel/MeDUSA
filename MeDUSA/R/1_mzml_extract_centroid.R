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
#'
#' @param files
#'   String: File or directory of mzmL files
#'   List: of mzML files (i.e. list.files(".",pattern = "qc_.*mzML"))
#' @param polarity
#'   0=negative, 1=positive, c(0,1)=both
#' @param params
#'   List: override any subset (or all) default parameters
#'   Default params Descriptions:
#'     "dbconn" connection object to a database. Or NULL to ignore database features
#'     "prebin_method" which sql_method(min, max, sum, avg) to combine intensities within the same mz bin per sample mzTs
#'     "postbin_method" which sql_method(min, max, sum, avg) to combine intensities within the same mz bin on the final mzObj
#'     "tolerance" tolerance for binning
#'     "timeSquash_method" which R_method(min, max, median, mean) to comibine intensities across scans
#'     "missingness_threshold" see mz_filter_missingness()
#'     "intensity_threashold" see mz_filter_low_intensity()
#'     "target_list" list of mzs to target
#'     "target_list_tolerance" tolerance of the mz target list
#'   Defaults params Values:
#'     "dbconn" = duckdb::dbConnect(duckdb::duckdb(paste0(local.output_dir(),'/',db_prefix,'_extract.duckdb'))),
#'     "prebin_method" = "max",
#'     "postbin_method" = "max",
#'     "tolerance" = 5e-6,
#'     "timeSquash_method" = mean,
#'     "missingness_threshold" = 1,
#'     "missingness_threshold" = .1,
#'     "intensity_threshold" = 1000,
#'     "target_list" = c(),
#'     "target_list_tolerance" = 5e-5
#' @examples
#'   mzml_extract_magic() : export all data in the current directory with default values
#'   mzml_extract_magic( getwd('data_dir'), polarity = 0, params = list("target_list" = c(100.111, 200.222, 300.333)))
#'      : export negative data, filtering on the targest provided, from the directory "data_dir".
#' @return MzObj & Database file in the output directory
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
      if (is.null(params$dbconn)) {
        stop("MeDUSA::mzml_extract_magic(): dbconn is currently required, do not set to NULL")
      }
      print(paste0("INFO: mzml_extract_magic : Extraction of [",pol_eng,"], beginning"))
      pbapply::pblapply(files, function(file) mzml_extract_file(file, polarity=pol, cl = NULL, magic=T, params=params, return_mzobj = F))
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
#'   parallel::threadCluster (optional).
#'      Note, this is untested after the introduction of duckDB, and the threads will likely compete in the DB. 
#'      It is suggested to choose threading, OR DuckDB (params = c( "dbconn"=NULL)). 
#' @param return_mzobj
#'   Boolean: Should this return mz_obj. True takes requires more memory, but is user friendly
#' @param params
#'   list of params, see mzml_extract_magic()
#' @examples
#'   mzml_extract_file("data_sample_001.mzml") : export data from "data_sample_001.mzml" with default parameters
#'   mzml_extract_file( "data_sample_001.mzml" , polarity = 0, params = list("target_list" = c(100.111, 200.222, 300.333)))
#'      : export negative data, filtering on the targest provided, from the file "data_sample_001.mzml".
#'   files <- list.files(path=getwd("data"), pattern="*.mzML");
#    mzT <- pbapply::pblapply(files, function(x) mzml_extract_file(x, polarity=0, cl = NULL, magic=F))
#'      : list mzml files in directory "data". Loop over mzml_extract, only negative values, single threaded, and without any additional filtering.
#'          This is good for troubleshooting, it will return a list of mzT objects with limited no additional processing.
#' @return MzT object
#' @export
mzml_extract_file <- function(file, polarity = "",  magic = T, cl = NULL, return_mzobj = T , params = NULL) {
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
    print(sprintf("INFO:MeDUSA::mzml_extract_file: Scans found in %s: %i",name, nrow(meta)))
    rts <-lapply(meta$retentionTime, function(x){c("mz",x)})
    l <- pbapply::pblapply(p, cl = cl, centroid.singleScan)
    for( i in seq_along(l) ){ 
      colnames(l[[i]]) <- rts[[i]] 
    }
    out <- extract.mz_format(l)
    if (magic) {
      print(sprintf("INFO:MeDUSA::mzml_extract_file: Binning & Filtering %s", name))
      out <- mzT_filtering( out,
                            prebin_method = params$prebin_method,
                            missingness_threshold = params$missingness_threshold,
                            intensity_threshold = params$intensity_threshold,
                            log_name = name)
      print(sprintf("INFO:MeDUSA::mzml_extract_file: TimeSquashing %s",name))
      out <- mzT_squash_time(out, timeSquash_method = params$timeSquash_method, output_name = name)
    }
    if (class(params$dbconn) == 'duckdb_connection'){
      duckdb::dbWriteTable(params$dbconn ,name , out,append = TRUE)
      if ( return_mzobj ){
        return(out)
      }
    } else {
      return(out)
    }
  } else {
    print(paste0("INFO:MeDUSA::mzml_extract_file: [", name, "] found in db, skipping"))
    if (return_mzobj) {
      out <- DBI::dbGetQuery(params$dbconn, paste0("select * from ", name))
    }
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
#' @param prebin_method \cr
#'   [R Method] : i.e. (mean, max, median)
#' @param missingness_threshold \cr
#'   Numeric   : If Decimal : Threshold % to blank data under. (as a decimal 10% = 0.1)
#'             : If Whole : Number of scans required to be nonzero ( Value required in at least 2 scans )
#' @param intensity_threshold \cr
#'   Numeric   : Lowest allowed intensity
#' @param log_name \cr
#'   String    : Identifier for log and plot outputs
#' @examples
#'   mzT_filtering(mzT) : Perform standard filtering
#'   mzT_filtering(mzT, prebin_method = 'mean', missingness_threshold = 3, intensity_threshold = 5000, log_name = "custom_filteing")
#'      : Filter using all custom values
#'   mzT_filtering(mzT, prebin_method = null, missingness_threashold = 1, intensity_threshold=0, log_name = "did_nothing")
#'      : Values required to skip all filtering
#' @return MzT object
#' @export
# Note missingness_threshold is very low
mzT_filtering <- function(mzT, prebin_method = 'max', missingness_threshold = 1, intensity_threshold = 1000, log_name = "" ){
  if ( !is.null(prebin_method)) {
    mzT <- mzT_binning(mzT, prebin_method,log_name)
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
#' @param timeSquash_method \cr
#'   [R method] : i.e. (mean, max, median)
#' @param ignore_zeros \cr
#'   Boolean: Should we set 0s to NA (to avoid effecting the math)
#' @param cl \cr
#'   parallel::threadCluster (optional).
#' @param output_name \cr
#'   String: Output name of the intensity column.
#' @examples
#'   mzT_squash_time(mzT) : Perform default time squashing
#'   mzT_squash_time(mzT, timeSquash_method = 'median', ignore_zeros = F)
#'      : Filter using custom values
#' @return MzObj of one sample column
#' @export
mzT_squash_time <- function(mzT, timeSquash_method = mean, ignore_zeros = T, cl = NULL, output_name = "intensity"){
  # This should be handled by filter low intesity
  if( ignore_zeros ){
    mzT[mzT == 0] <- NA
  }
  squash <- pbapply::pbapply(dplyr::select(mzT, -mz),1, timeSquash_method, na.rm=T, cl=cl)
  squash[is.na(squash)] <- 0
  out <- data.table::as.data.table(cbind(mzT$mz, squash))
  names(out) <- c("mz", output_name)
  out
}

# ***  -----------------------------------------------------
#' Combine mzs within a tolerance of other mzs
#'
#' @param mzT \cr
#'   MzT : MzT object (also will work on MzObj)
#' @param method \cr
#'   [Sql method] : i.e. (avg, max, median)
#' @param tolerance \cr
#'   Tolerance for binning
#' @param log_name \cr
#'   String    : Identifier for db, log and plot outputs
#' @examples
#'   mzT_binning(mzT) : Perform default binning
#'   mzT_binning(mzT, method = 'mean', log_name = "custom_values", tolerance = 1e-8)
#'      : Binning using custom values
#' @return MzObj of combined mz rows
#' @export
mzT_binning <- function(mzT, method = 'max', log_name = "", tolernace = 5e-6){
  db_name <- paste0("Binning_", log_name)
  db_file <- paste0(local.output_dir(),'/tmp-',db_name,'.duckdb')
  bin <- binning(mzT$mz, tolerance = 5e-6)
  local.mz_log_removed_rows((bin), unique(bin), db_name)
  mzT$mz <- bin
  rm(bin); gc()
  #TODO see if we can eliminate this check (scanTime as string from extract_file gets weird)
  if (regexpr("extract_magic.*",log_name)){
    cols <- colnames(dplyr::select(mzT, -mz))
  } else {
    cols <- as.numeric(colnames(dplyr::select(mzT, -mz)))
  }
  cols <- sprintf('%s("%s")', method, paste(cols ,collapse=paste0('"),',method,'("')))
  #TODO: investigate possible duplications
  query <- paste0("CREATE TABLE final_tmp AS (SELECT mz, ",cols," FROM ", db_name, " GROUP BY mz ORDER by mz)" )
  #Create tmp db to perform binning much quicker via SQL (and with limited memory usage)
  dbconn <- DBI::dbConnect(duckdb::duckdb(db_file))
  duckdb::duckdb_register(dbconn, db_name, mzT, overwrite = T)
  rm(mzT); gc()
  DBI::dbExecute(dbconn, query)
  mzT <- DBI::dbGetQuery(dbconn, "SELECT * FROM final_tmp")
  DBI::dbExecute(dbconn, "DROP TABLE final_tmp")
  DBI::dbDisconnect(dbconn)
  unlink(db_file)
  colnames(mzT)<- gsub(paste0(method,"\\((.*)\\)"),"\\1",colnames(mzT))
  mzT[mzT < 0] <- 0
  mzT[is.na(mzT)] <- 0
  mzT
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
    warning("ERROR:MeDUSA::mzml_extract_magic: Cannot find mzML files in given files.")
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
    warning("WARN:MeDUSA::fill_defaults: Binning method should be SQL aggregate function (max, min, sum, avg)")
  }

  out
}