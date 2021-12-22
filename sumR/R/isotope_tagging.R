#-----------------------------------------------------------------
#functions to use in Isotope_tagging
#' @title Dalton error calculation
#' @description This function computes the error in Dalton
#' @usage ppm_to_dalton(mass)
#' @param mass Mass of interest in m/z
#' @param ppm  Parts Per Million Tolerance
ppm_to_dalton <- function(mass, ppm) {
  da_error  <- (mass*ppm) / 1e6
  return(da_error)}

#' @title Probability mono-isotopic peak calculation for carbon elements
#' @description This function computes the probability of the mono-isotopic peak
#' @usage Prob_M(n)
#' @param n number of atoms
Prob_M <- function(n) {
  prob_m  <- (1 - 0.0107)^n
  return(prob_m)}

#' @title Probability isotopic peak calculation for carbon elements
#' @description This function computes probability isotope calculation
#' @usage Prob_M_1(n)
#' @param n number of atoms
Prob_M_1 <- function(n) {
  prob_m_plus1  <- n * (0.0107 / (1 - 0.0107)) * Prob_M(n)
  return(prob_m_plus1)}

#' @title Estimation of number of carbon in Alkyne
#' @description Calculates the number of carbons in C(n)H(2n-2)
#' @usage Alkyne_n(n)
#' @param mz Mass of interest in m/z
Alkyne_n <- function(mz) {
  n  <- round((mz - 2) / 14)
  return(n)}

#' @title Intensity ratio for carbon
#' @description Calculates the intensity ratio for carbon C(n)H(2n-2) within a specified range
#' @usage Alkyne_n(n)
#' @param mz Mass of interest in m/z
#' @param by_step "by" in sequence range
Intensity_ratio_range <- function(mz,by_step){
  seq_vec  <- seq(from = min(mz) , to = max(mz), by  = by_step)
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

#' @title Filtering for intensity ratio
#' @description Filtering for intensity ratio within each range acccording to Intensity ratio vector
#' @usage subset_matrix(df, mz , Intensity_ratio_vec, by_step)
#' @param df Data frame
#' @param mz Mass
#' @param Intensity_ratio_vec Intensity ratio vector for carbon
#' @param by_step "by" in sequence range
#' @importFrom dplyr select mutate case_when
subset_matrix_filter <- function(df, mz , Intensity_ratio_vec, by_step ){
  seq_vec  <- seq(from = min(mz) , to = max(mz), by  = by_step)
  result <- c()
  for (j in 1:length(Intensity_ratio_vec)) {
    for (i in 1:length(seq_vec)) {
      if (i == j) {
        df <- df %>%
          select(mz, intensity) %>%
          mutate(
            temp_intensity = case_when(
              mz >= seq_vec[i] & mz < seq_vec[i + 1] ~ intensity ))

        iso_ratio_mat <- outer(df$temp_intensity , df$temp_intensity , `/`)
        index_ratio_mat <- which(iso_ratio_mat > Intensity_ratio_vec[j] & iso_ratio_mat < Intensity_ratio_vec[j + 1] , arr.ind = TRUE )
        result[[i]] <- index_ratio_mat
      }
    }
  }
  return(result)
}


#' @title Mass spectrometer visualization
#' @description Gives a graphical representation of the mono-isotopic and isotopic peaks
#' @usage mass_spec_plot(final_df_vis)
#' @param final_df_vis dataframe of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  plotly plot_ly highlight_key add_markers add_segments  layout rangeslider ggplotly highlight
mass_spec_plot <- function(final_df_vis){
  d <- highlight_key(final_df_vis, ~id )

  fig <- plot_ly(data = d, x = ~mz, y = ~intensity, color = ~isotopic_status,colors = "Set1",
                 text = ~paste("<b>m/z </b>: ",mz, '<br><b>id</b>:', id), hoverinfo = "text")
  fig <- fig %>%
    add_markers() %>%
    add_segments(
      data = final_df_vis, x = ~mz, xend = ~mz,
      y = 0, yend = ~intensity, color = I("grey"), showlegend = FALSE)  %>%
    rangeslider(final_df_vis$mz[1], max(final_df_vis$mz))

  gg <- ggplotly(fig, tooltip = "id")
  gg <- highlight(gg, dynamic = TRUE)
  return(gg)
}



#---------------------------------------------------------------------------------------------------
#' @title Isotope tagging
#' @description Tags isotopes based on intensity ratio and difference in mass
#' @usage Isotope_tagging(df, iso_diff , ppm , by_step )
#' @param df Dataframe must contain columns c("mz", "intensity")
#' @param iso_diff_da Mass difference between mol ion and isotope
#' @param ppm Parts Per Million Tolerance
#' @param by_step "by" in sequence mz range
#' @importFrom  data.table as.data.table fintersect
#' @importFrom  dplyr tibble filter select mutate
isotope_tagging <- function(df , iso_diff_da, ppm , by_step ){
#-----------------------------------------------------------------------------------
## Start body part of function BAIT
DT <- as.data.table(df)
#setting mz & intensity to vectors
mz_vector <-  DT[, mz]
int_vector <- DT[, intensity]

#creating new df with mass_error , lower bound and upperbound
  df <- mutate(df,
               mass_error = ppm_to_dalton(mz, ppm),
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
  Intensity_ratio_vec <- Intensity_ratio_range(mz_vector, by_step )
  list_index_ratio_matrix <- subset_matrix_filter(df,mz_vector, Intensity_ratio_vec, by_step)
  index_mat_ratio <- do.call(rbind, list_index_ratio_matrix)

  #transforming index matrix to datatable
  index_dt_min <- as.data.table(index_min_mz )
  index_dt_max <- as.data.table(index_max_mz)
  index_dt_ratio <- as.data.table(index_mat_ratio)


  #combining innerjoin , has to fulfill (1)isotope_da +- tolerance "Da error"
  #                                     (2)iso-ratio <= than specific ratio within mz range
  index_dt_range <- fintersect(index_dt_min, index_dt_max , all = TRUE)
  index_matches <- fintersect(index_dt_range, index_dt_ratio , all = TRUE)


  #creating dataframe
  iso_df <- tibble(mol_ion = mz_vector[index_matches[,col]],
                   isotope = mz_vector[index_matches[,row]],
                   pair = paste("[",mol_ion,",",isotope,"]"),
                   diff = isotope - mol_ion)
  iso_df$id <- paste0('id-', 1:nrow(iso_df))

  #tagging
  final_df <- select(mutate(df,
                            isotopic_status = case_when(
                              mz_vector %in% iso_df$mol_ion ~ "molecular_ion",
                              mz_vector %in% iso_df$isotope ~ "c13_isotope")),
                     isotopic_status ,
                     everything())


  final_df_iso <- merge(final_df, iso_df ,by.x = "mz", by.y = c("isotope"))
  final_df_mon <- merge(final_df, iso_df ,by.x = "mz", by.y = c("mol_ion"))

  final_df_iso$mol_ion <- NULL
  final_df_mon$isotope <- NULL
  final_df_vis <- rbind(final_df_mon,final_df_iso)


  final_df <- merge(final_df, final_df_vis, by = "mz", all = TRUE)
  final_df <- final_df %>%
    rename(
      Isotopic_Status = isotopic_status.x,
      Difference_Da = diff,
      Mol_Isotope_Pair = pair) %>%
    mutate(Confidence_ppm = paste("(",round(min_error.x,2),",",round(max_error.x,2),")"))

  final_df <- final_df[, c(2,1,15,12,13,14)]
  final_df

  #plotting
  print(mass_spec_plot(final_df_vis))

  return(final_df)
}





