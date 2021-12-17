#functions to use in Isotope_tagging
#' @title Dalton error calculation
#' @description This function computes the error in Dalton
#' @usage ppm_to_dalton(mass)
#' @param mass Mass of interest in m/z
#' @param ppm  Parts Per Million Tolerance
#' @example
#'
ppm_to_dalton <- function(mass) {
  da_error  <- (mass*ppm) / 1e6
  return(da_error)}


#' @title Probability monoisotopic calculation carbon
#' @description This function probability of the monoisotopic
#' @usage Prob_M(n)
#' @param n number of atoms
#' @example
#'
Prob_M <- function(n) {
  prob_m  <- (1 - 0.0107)^n
  return(prob_m)}

#' @title Probability isotope calculation  carbon
#' @description This function computes Probability isotope calculation
#' @usage Prob_M_1(n)
#' @param n number of atoms
#' @example
#'
Prob_M_1 <- function(n) {
  prob_m_plus1  <- n * (0.0107 / (1 - 0.0107)) * Prob_M(n)
  return(prob_m_plus1)}

#' @title Number of carbon in Alkyne
#' @description calculates the number of carbons in C(n)H(2n-2)
#' @usage Alkyne_n(n)
#' @param mz Mass of interest in m/z
#' @example
#'
Alkyne_n <- function(mz) {
  n  <- round((mz - 2) / 14)
  return(n)}

#' @title Intensity ratio for carbon
#' @description calculates the number of carbons in C(n)H(2n-2)
#' @usage Alkyne_n(n)
#' @param mz Mass of interest in m/z
#' @example
#'
Intensity_ratio <- function(mz, by_step) {
  seq_vec  <- seq(from = mz[1], to = length(mz), by = by_step)
  for (i in 1:length(seq_vec)) {
    n <- Alkyne_n(seq_list)
    Prob_M_res <- Prob_M(n)
    Prob_M_1_res <- Prob_M_1(n)
    Intensity_ratio <- (Prob_M_1_res / Prob_M_res)
  }
  return(Intensity_ratio)}


#' @title Filtering for intensity ratio
#' @description Filtering for intensity ratio within each range
#' @usage subset_matrix(mat, n_rows, n_cols)
#' @param mat matrix with isotopic ratio
#' @param n_rows number of rows for submatrix
#' @param n_cols number of columns for submatrix
#' @example
#'
subset_matrix <- function(mat, n_rows, n_cols){
  # Loop over the matrix, increase the number of rows with n_rows
  for (row in seq(from = 1, to = nrow(mat), by = n_rows)) {
    # Determine maximum row-index
    max_row <- row + n_rows - 1

    # Determine maximum column-index
    max_col <- row + n_cols - 1


    # Check that max row is not bigger than
    # the number of rows in the matrix
    if (max_row <= nrow(mat)) {

      # Subset the matrix
      sub_matrix <- mat[row:max_row, row:max_col]
      # Add any futher modification
      sub_matrix_list <-  c(sub_matrix_list,sub_matrix)
      return(sub_matrix_list)
    }

  }
}

#' @title Dalton error calculation
#' @description This function computes the error in Dalton
#' @usage ppm_to_dalton(mass)
#' @param mass Mass of interest in m/z
#' @param ppm  Parts Per Million Tolerance
#' @example
#'
ppm_to_dalton <- function(mass) {
  da_error  <- (mass*ppm) / 1e6
  return(da_error)}









#---------------------------------------------------------------------------------------------------



#' @title Isotope tagging
#' @description Tags isotopes based on intensity ratio and
#' @usage
#' @param df Dataframe
#' @param iso_diff_da Mass difference between mol ion and isotope
#' @param iso_ratio   Ratio between I1/I2
#' @param ppm Parts Per Million Tolerance
#' @importFrom  data.table as.data.table finstersect
#' @importFrom  dplyr tibble filter select mutate
isotope_tagging <- function(df, iso_diff_da, iso_ratio, ppm){

  #-----------------------------------------------------------------------------------
  ## Start body part of function BAIT
  DT <- as.data.table(df)
  #setting mz & intensity to vectors
  mz_vector <-  DT[, mz]
  int_vector <- DT[, intensity]

  #creating new df with mass_error , lower bound and upperbound
  df <- mutate(df,
               mass_error = ppm_to_dalton(mz),
               min_error  = mz - mass_error,
               max_error  = mz + mass_error)

  #assigning minimum error and max error
  min_error_vec <- df$min_error
  max_error_vec <- df$max_error

  #computing maximum difference and minimum difference
  max_diff_mat <- abs(outer(max_error_vec, min_error_vec, `-`))
  min_diff_mat <- abs(outer(min_error_vec, max_error_vec, `-`))

  #filtering + indexing max and min difference
  index_max_mz <- which(max_diff_mat >= isotope_da,
                        arr.ind = TRUE )
  index_min_mz <- which(min_diff_mat <= isotope_da & min_diff_mat > 0,
                        arr.ind = TRUE )


  #-------------------------------------------------------------------------------
  #selecting every isotope candidate with ratio <0.5
  iso_ratio_mat <- outer(int_vector , int_vector , `/`)
  index_mat_ratio <- which(iso_ratio_mat <= iso_ratio, arr.ind = TRUE)






  #transforming index matrix to datatable
  index_dt_min <- as.data.table(index_min_mz )
  index_dt_max <- as.data.table(index_max_mz)
  index_dt_ratio <- as.data.table(index_mat_ratio)


  #combining innerjoin , has to fulfill (1)isotope_da +- tolerance "Da error"
  #                                     (2)iso-ratio <0.5
  index_dt_range <- fintersect(index_dt_min, index_dt_max , all = TRUE)
  index_matches <- fintersect(index_dt_range, index_dt_ratio , all = TRUE)


  #creating column  => will need to be adjusted at the end
  iso_df <- tibble(mol_ion = mz_vector[index_matches[,col]],
                   isotope = mz_vector[index_matches[,row]],
                   pair = paste("pair(",mol_ion,",",isotope,")"),
                   diff = isotope - mol_ion)
  #filtering
  iso_df_f <- filter(iso_df, diff > 0)

  #tagging
  final_df <- select(mutate(df,
                            isotopic_status = case_when(
                              mz_vector %in% iso_df_f$mol_ion ~ "molecular_ion",
                              mz_vector %in% iso_df_f$isotope ~ "c13_isotope")),
                     isotopic_status ,
                     everything() )

  final_df <- as.data.frame(final_df)

  return(final_df)  }
