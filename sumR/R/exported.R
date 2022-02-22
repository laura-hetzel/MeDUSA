#' Iso hippo
#' 
#' @description iso_hippo is a function that takes 1 argument, and 2 optional ones
#' Similarly to kopi luwak, iso hippo is a retired, senile, hippopotamus
#' that grazes on huge chunks of data and poops isotopes and mol.ions
#'
#' Ahmed's notes: Should expand on it and allow it to work with other isotopes +
#' merge it with the repeating units function, also have it check intensity values
#' #to make sure if its actually an isotope of the parent molecule or not
#' @param df is the dataframe on which hippo eats, it has to have an mz column
#' @param isotope_da is an optional argument, specifying the differences
#' between mol ions and isotopes in Da
#' @param iso_tolerance optional argument: tolerance in Da for isotope differences
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
#' @examples
#'
#'
#' @export
iso_hippo <- function(df, isotope_da = 1.0034, iso_tolerance = 0.0034) {
  mz_vector <- as.numeric(df$mz) #incase mzs are strings
  ooutput <- abs(outer(mz_vector, mz_vector, `-`))
  index_matrix <- which(ooutput <= isotope_da + iso_tolerance &
                          ooutput >= isotope_da - iso_tolerance,
                        arr.ind = TRUE)
  iso_df <- tibble(mol_ion = mz_vector[index_matrix[, 1]],
                   isotope = mz_vector[index_matrix[, 2]],
                   diff = mol_ion - isotope)
  iso_df_f <- filter(iso_df, diff > 0)
  final_df <- select(mutate(df, isotopic_status =
                              case_when(mz_vector %in% iso_df_f$mol_ion ~ "molecular_ion",
                                        mz_vector %in% iso_df_f$isotope ~ "c13_isotope")),
                     isotopic_status, everything())
  return(final_df)
}


#' Title
#'
#' @param rnamed_cnamed_df data.frame with rownamed mz,colnamed samples
#' @param file_path desire output path, should be a string
#' @param extension desired custom extension to files (aside from .txt), string
#' 
#' @description takes the df, separates it into mz,int files for each sample (1st loop)
#' export each file with desired file path,and custom extension (2nd loop)
#'
#' @return
#' @export
#' @importFrom utils write.table
#'
#' @examples
FUN_text_exporter <- function(rnamed_cnamed_df, file_path, extension) {
  all_spectra <- list()
  c_names <- colnames(rnamed_cnamed_df)
  r_names <- row.names(rnamed_cnamed_df)
  message("column /row names extracted , commencing 1st boop de loop")
  for (i in seq_along(c_names)) {
    assign(paste0(c_names[i],'_BG'),
           as.data.frame(cbind(mz = as.numeric(r_names),
                               int = rnamed_cnamed_df[ , c_names[i]])))
  }
  message("first loop successiful, commecning 2nd boop de loop")

  all_spectra <- lapply(ls()[grepl(pattern = '_BG', ls())], get, envir = environment())
  message(paste("Length of file list: ",length(all_spectra)))

  filenames <- c_names
  for (i in seq_along(filenames)) {
    write.table(subset(all_spectra[[i]], all_spectra[[i]][2] > 0),
                file = paste0(file_path,filenames[i],
                              extension,'.txt'),
                col.names = FALSE, row.names = FALSE, sep = '\t')
  }
  message("everything is done")
}

