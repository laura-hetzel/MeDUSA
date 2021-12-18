#-----------------------------------------------------------------
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
#' @description This function computes the probability of the monoisotopic atom
#' @usage Prob_M(n)
#' @param n number of carbon atoms
#' @example
#'
Prob_M <- function(n) {
  prob_m  <- (1 - 0.0107)^n
  return(prob_m)}

#' @title Probability isotope calculation  carbon
#' @description This function computes the probability of the isotope atom
#' @usage Prob_M_1(n)
#' @param n number of carbon atoms
#' @example
#'
Prob_M_1 <- function(n) {
  prob_m_plus1  <- n * (0.0107 / (1 - 0.0107)) * Prob_M(n)
  return(prob_m_plus1)}

#' @title Number of carbons in Alkyne
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
#' @usage Intensity_ratio_range(df, by_step)
#' @param df dataframe that contains mz and intensity
#' @param mz range for evaluation
#' @example
#'
Intensity_ratio_range <- function(df,by_step){
  seq_vec  <- seq(from = min(df$mz) , to = max(df$mz), by  = by_step)
  Intensity_ratio_vec = numeric(length(seq_vec) )
  for (i in 1:length(seq_vec)) {
    n <- Alkyne_n(seq_vec[i])
    Prob_M_res <- Prob_M(n)
    Prob_M_1_res <- Prob_M_1(n)
    Intensity_ratio <- (Prob_M_1_res / Prob_M_res)
    Intensity_ratio_vec[i] <- Intensity_ratio
  }
  return(Intensity_ratio_vec)
}



#' @title Filtering for intensity ratio within mz range
#' @description Filtering for intensity ratio within each range of mz
#' @usage subset_matrix_filter(df , Intensity_ratio_vec, by_step )
#' @param df dataframe that contains mz and intensity
#' @param Intensity_ratio_vec vector of intensity ratio's for each mz range
#' @param by_step mz range for evaluation
#' @importFrom  dplyr select mutate case_when
#' @example
#'
subset_matrix_filter <- function(df , Intensity_ratio_vec, by_step ){
  seq_vec  <- seq(from = min(df$mz) , to = max(df$mz), by  = by_step)
  result <- c()
  for (j in 1:length(Intensity_ratio_vec)) {
    for (i in 1:length(seq_vec)) {
      if (i == j) {
        df <- df %>%
          select(mz, intensity) %>%
          mutate(
            temp_intensity = case_when(
              mz >= seq_vec[i] & mz < seq_vec[i + 1] ~ intensity
            )
          )
        iso_ratio_mat <- outer(df$temp_intensity , df$temp_intensity , `/`)
        index_ratio_mat <- which(iso_ratio_mat > Intensity_ratio_vec[j] & iso_ratio_mat < Intensity_ratio_vec[j + 1] , arr.ind = TRUE )
        result[[i]] <- index_ratio_mat
      }
    }
  }
  return(result)
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
#' @param by_step The mz range for evaluation
#' @importFrom  data.table as.data.table finstersect
#' @importFrom  dplyr tibble filter select mutate case_when
#'
isotope_tagging <- function(df, by_step, ppm, iso_diff_da ){
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
  index_max_mz <- which(max_diff_mat >= iso_diff_da,
                        arr.ind = TRUE )
  index_min_mz <- which(min_diff_mat <= iso_diff_da & min_diff_mat > 0,
                        arr.ind = TRUE )


  #-------------------------------------------------------------------------------
  #selecting every isotope candidate with ratio iso-ratio <= than specific ratio within mz range
  Intensity_ratio_vec <- Intensity_ratio_range(df$mz, by_step )
  list_index_ratio_matrix <- subset_matrix_filter(df, Intensity_ratio_vec, by_step)
  index_mat_ratio <- do.call(rbind, list_index_ratio_matrix)

  #transforming index matrix to datatable
  index_dt_min <- as.data.table(index_min_mz )
  index_dt_max <- as.data.table(index_max_mz)
  index_dt_ratio <- as.data.table(index_mat_ratio)


  #combining innerjoin , has to fulfill (1)isotope_da +- tolerance "Da error"
  #                                     (2)iso-ratio <= than specific ratio within mz range
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
  return(final_df)
}



