#' @title Perform a function without message or output
#' @description Use this function to completely block any output messages
#' a function has.
#' @details Some functions use either a form of logging, messaging or `cat` to
#' produce some output. However, the output is not always necessary or doesn't
#' fit the output desired in sumR. Enclosing the function with a
#' `quiet(myfunc)` will ignore any output `myfunc` will produce.
#' @returns When used as `quiet(myfunc)` this function will return the output
#' of `myfunc`.
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
#' @details This function (together with `saveExperiment`) aids in
#' loading/storing experiments at different post-processing stages. It can also
#' be shared across users.
#' @returns A SummarizedExperiment that was saved during a previous run.
#' @param path Character path to the location where the RDS file will be stored.
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

#' @title Save a sumR SummarizedExperiment to an RDS file
#' @description This is a helper function to store the SummarizedExperiment
#' as an RDS file after validation.
#' @details Saving an experiment can be useful for multiple purposes,
#' including making comparisons and sharing. This function performs a few
#' validation checks before it creates the directory of `path` if it does not
#' exist yet. After creation, the [saveRDS] function is used to store a
#' compressed file of the SummarizedExperiment to the disk.
#'
#' Similarly, loading an experiment can be done with the [loadExperiment]
#' function.
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

#' @title Combine SummarizedExperiments by row(s)
#' @description This function intersects 2 or more SummarizedExperiments into
#' a single object. This can be used to combine positive and negative polarity
#' into one SummarizedExperiment.
#' @details While in an experiment, often both polarities are measured, each
#' polarity must be analyzed separately. This ensures that peaks aren't aligned
#' across polarities. However, before modelling, the polarities should be
#' combined in order to find the combination of peaks that are able to
#' separate phenotypes. This function can bind two polarities that use the
#' same samples together. This acts similarly as [base::rbind], by binding the
#' [rowData] and [assay]s of the [SummarizedExperiment].
#' @returns A SummarizedExperiment if the experiments supplied are valid.
#' Otherwise the list of objects is returned.
#' @param exps A list of SummarizedExperiments to be combined
#' @export
combineExperiments <- function(exps){
  if (!all(sapply(exps, validateExperiment))) return(exps)

  idx <- Reduce(intersect, lapply(exps, colnames), init = colnames(exps[[1]]))
  exp <- do.call(rbind, lapply(exps, function(x) x[, idx]))
  rownames(exp) <- 1:nrow(exp)
  exp
}

#' @title Convert vendor files to mzML using Proteowizard
#' @description Vendor files obtained with mass spectrometry have their own
#' (compressed) file format. This function helps in converting these files
#' to the open mzML format. This function requires Proteowizards __msconvert__
#' has been installed on your system.
#' @details mzML is the most common format for mass spectrometry files. Here,
#' using the help of Proteowizards msconvert, vendor-formatted files can be
#' converted to mzML. By default, files are stored in the same folder as the
#' vendor-formatted files. Additionally, the option `--simAsSpectra` is used
#' to convert Selected Ion Chromatograms as spectra.
#' @returns A character vector of path(s) to the stored mzML files. These
#' are often used as input to [extractPeaks]
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
