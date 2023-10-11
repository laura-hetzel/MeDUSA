# *** IsotopeTagging -----------------------------------------------------
#' MZ-OBJ IsotopeTagging
#'
#' Hunts for isotopes. Defaults to Carbon13
#'   You may use the return to compare intensities of isotopes.
#'   input_mz_obj[input_mz_obj$mz %in% output[['100.12345']],]
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param iso_target \cr
#'   Float : Dalton diff of isotope to hunt
#' @param iso_iter \cr
#'   Integer     : How many 'iterations' of iso_target should we look
#' @param tol_ppm \cr
#'   Float   : Overwrite csv tolerance, with single value ppm
#' @param cores \cr
#'   Integer : How many threads could a multi-threader thread, if a multi-threader could thread threads?
#'
#' Dependencies : pbapply, dplyr
#' @return Returns a List of MZ that are isotopes of each other
#' @export

mz_isotope_hunter <- function(input_mz_obj, iso_target = 1.0034, iso_iter = 5, tol_ppm = 5e-6, cores = 4, ...){
  input_mz_obj$tol <- input_mz_obj$mz - (input_mz_obj$mz/( tol_ppm + 1 ))
  sub_set <- dplyr::select(input_mz_obj, c(mz, tol))
  out <- list()

  cl <- local.export_thread_env(cores, deparse(sys.calls()[[sys.nframe()]]))
  tryCatch({
    out <- pbapply::pbapply(sub_set, 1,cl=cl,function(mz_target){
      hunt_zone  <- dplyr::filter(sub_set, (mz_target[1] - iso_target * (iso_iter + 0.1)) < mz
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
    if (cores > 1 || !is.null(cl)) {
      parallel::stopCluster(cl)
      showConnections()
    }
  })
}

#mz_isotope_flatten <- function(input_mz_obj, isotope_list, method = max, cores = 4){
#  pblapply(input_mz_obj, function(isotope_list, x, method = method){
#    if x$mz %in% names(isotope_list){
#      x <- apply(select(input_mz_obj[input_mz_obj$mz %in% x,], -mz),2, method)
#  })
#}


#TODO FLATTEN?
