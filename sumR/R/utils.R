#' @title Perform a function without message or output
#' @param x function to be executed
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' @title Load a sumR Experiment from a file
#' @description SummarizedExperiments can be saved during any step after
#' alignment using the `saveExperiment` function. This function reads the
#' resulting RDS file that saveExperiment produces. This allows for saving
#' work and continuing it later.
#' @param path Character path to the location where the RDS file is stored.
#' @export
loadExperiment <- function(path){
  tryCatch({
    exp <- readRDS(path)
    if (validateExperiment(exp)) {
      return(exp)
    }
    warning("Invalid Experiment, returning NULL")
    return(NULL)
  }, error = function(x){
    message("Error while trying to read the RDS file")
  })
}

#' @title Save a sumR Experiment to a file
#' @param exp SummarizedExperiment obtained after alignment
#' @param path Character of the location where to store the Experiment as RDS
#' file.
#' @export
saveExperiment <- function(exp, path){
  path <- file.path(path)
  if (validateExperiment(exp)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(exp, file = path)
  } else {
    warning("Invalid Experiment, aborted saving")
  }
}

#' @title Combine SummarizedExperiments by row
#' @param exps A list of SummarizedExperiments to be combined
#' @export
combineExperiments <- function(exps){
  if (!all(sapply(exps, validateExperiment))) return(NULL)

  idx <- Reduce(intersect, lapply(exps, colnames), init = colnames(exps[[1]]))
  exp <- do.call(rbind, lapply(exps, function(x) x[, idx]))
  rownames(exp) <- 1:nrow(exp)
  exp
}

#' @title Convert vendor files to mzML using Proteowizard
#' @param folder Character path to the folder containing vendor-formatted files,
#' e.g. wiff from Sciex.
#' @param output Folder where to store the resulting mzML files.
#' @param rt A 2-sized vector of rt-min and rt-max forming a region to extract
#' from the vendor-formatted file. Can be used to extract a region of interest.
#' @param options Additional options for msConvert. Defaults to
#' `--simAsSpectra`. Any options given will overwrite this option, unless it is
#' added to the options.
#' @export
rawToMzml <- function(folder, output = getwd(), rt = NULL,
                      options = "--simAsSpectra"){

  options <- paste(options, collapse = " ")
  files <- list.files(file.path(folder), full.names = T)
  files <- files[!grepl(".mzML$", files)]
  if (length(files) == 0) stop("No files found to convert")

  output <- file.path(output)
  if (!dir.exists(output)) dir.create(output, recursive = T, showWarnings = F)

  key <- "Directory\\shell\\Open with MSConvertGUI\\command"
  reg <- tryCatch(utils::readRegistry(key, "HCR"), error = function(e) NULL)
  msconvert <- file.path(dirname(sub("\"([^\"]*)\".*", "\\1", reg[[1]])), "msconvert.exe")
  if (is.null(msconvert)) return(NULL)

  rt <- ifelse(is.null(rt), "", sprintf('--filter "scanTime [%s, %s]"', rt[1], rt[2]))

  pbapply::pblapply(files, function(file){
    command <- sprintf('"%s" "%s" %s %s -o "%s"',
                       file.path(msconvert), file, rt, options, output)
    system(command, show.output.on.console = F)
  })
  list.files(output, full.names = T, pattern = ".mzML")
}
