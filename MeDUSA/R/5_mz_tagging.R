# *** IsotopeTagging -----------------------------------------------------
#' Find Isotopes (default c13)
#'
#' @description
#' In mass spectrometry it is well known that every m/z value measured may not
#' be a unique feature, largely due to isotopes. The mz_tag_isotope_hunter
#' function analyzes all m/z features and identifies values that are likely an
#' isotope instead of a unique feature. It is up to the user to remove these
#' m/z values if desired.
#'
#' MZ-OBJ IsotopeTagging
#'
#' Hunts for isotopes. Defaults to Carbon13
#'   You may use the return to compare intensities of isotopes; mind trunctation.
#'   input_mz_obj[input_mz_obj$mz %in% output[['100.12345']],]
#'
#' @param input_mz \cr
#'   DataFrame : Input MZ-Obj (Can take a full mz_obj or mm_obj$mz)
#' @param iso_target \cr
#'   Float : Mass of isotope to hunt.
#' @param iso_iter \cr
#'   Integer     : How many 'iterations' of iso_target should we look
#' @param tol_ppm \cr
#'   Float   : Overwrite csv tolerance, with single value ppm
#' @param cores \cr
#'   Integer : How many threads could a multi-threader thread, if a multi-threader could thread threads?
#'
#' Dependencies : pbapply, dplyr
#'
#' @examples
#'   mz_tag_isotope_hunter(input_mz) : Find isotopes of carbon 13
#'   mz_tag_isotope_hunter(input_mz, iso_target = 100.123, iso_iter = 8, tol_ppm = 1e-7)
#'     : Find isotopes on a different isotope, with a wider range, and stricter tolerance
#'
#' @return Returns a List of MZ that are isotopes of each other
#' @export

mz_tag_isotope_hunter <- function(input_mz, iso_target = 1.0034, iso_iter = 5, tol_ppm = 5e-6, cores = 4, ...){
  if( !is.null(ncol(input_mz)) ){
    input_mz <- input_mz$mz
  }
  input_mz <- data.frame(mz  = input_mz,
                         tol = input_mz - (input_mz/( tol_ppm + 1 )))
  out <- list()

  cl <- local.export_thread_env(cores, environment())
  tryCatch({
    out <- pbapply::pbapply(input_mz, 1,cl=cl,function(mz_target){
      hunt_zone  <- dplyr::filter(input_mz, (mz_target[1] - iso_target * (iso_iter + 0.1)) < mz
                            & mz <  (mz_target[1] + iso_target * (iso_iter + 0.1)))
      x <- hunt_zone$mz[(((abs(hunt_zone$mz - mz_target[1]) + mz_target[2]) %% iso_target )
                            - mz_target[2] ) < mz_target[2]]
      out[[mz_target[1]]] <- x
    })

  out <- unique(out[ lapply(out, length) > 1 ])
  names(out) <- lapply(out, `[[`, 1)

  return(out)
  },
  finally={
    local.kill_threads(cl)
  })
}

#mz_isotope_flatten <- function(input_mz_obj, isotope_list, method = max, cores = 4){
#  pblapply(input_mz_obj, function(isotope_list, x, method = method){
#    if x$mz %in% names(isotope_list){
#      x <- apply(select(input_mz_obj[input_mz_obj$mz %in% x,], -mz),2, method)
#  })
#}

#TODO FLATTEN?
