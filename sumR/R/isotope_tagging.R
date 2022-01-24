#-----------------------------------------------------------------
#functions to use in Isotope_tagging
#' @title Dalton error calculation
#' @description This function computes the error in Dalton
#' @param mass Mass of interest in m/z
#' @param ppm  Parts Per Million Tolerance
#' @usage ppm_to_dalton(mass, ppm)
ppm_to_dalton <- function(mass, ppm) {
  da_error  <- (mass*ppm) / 1e6
  return(da_error)}



#' @title Adding intensities
#' @description A list of annotation information of isotopes
#' @param df input dataframe
#' @param iso_df data frame adducts
#' @usage add_intensities(df, iso_df)
#' @importFrom dplyr filter 
add_intensities <- function(df, iso_df){
  add_mol <- merge(filter(df, df$mz %in% iso_df$mol_ion),iso_df, by.x = "mz", by.y = "mol_ion", all.x = TRUE)
  add_mol <- add_mol[,c(1,2)]
  add_mol$intensity_mol <- add_mol$intensity
  add_mol$mz_mol <- add_mol$mz
  add_mol$intensity <- NULL
  add_mol$mz <- NULL
  
  
  add_iso <- merge(filter(df, df$mz %in% iso_df$isotope),iso_df, by.x = "mz", by.y = "isotope", all.x = TRUE)
  add_iso <- add_iso[,c(1,2)]
  add_iso$intensity_iso <- add_iso$intensity
  add_iso$mz_iso <- add_iso$mz
  add_iso$intensity <- NULL
  add_iso$mz <-NULL
  
  iso_df_fin <- cbind(add_mol,add_iso)
  
  
  return(iso_df_fin)
}

#' @title Annotation list for isotopes
#' @description A list of annotation information of isotopes
#' @param iso_df_fin dataframe filtered for isotopes with possible isotopes
#' @param ppm error part per million
#' @param  z charge values
#' @usage isotope_molecules(iso_df_fin, ppm , z )
#' @importFrom Rdisop decomposeIsotopes
isotope_molecules <- function(iso_df_fin, ppm ,z){
  isotope_molecules_list <- NULL
  length_iso <- nrow(iso_df_fin)
  for (i in 1:length_iso) {
    masses <- c(iso_df_fin$mz_mol[i], iso_df_fin$mz_iso[i])
    intensities <- c(100, ((iso_df_fin$intensity_iso[i]/iso_df_fin$intensity_mol[i])*100) )
    isotope_molecules_list[[i]] <- decomposeIsotopes(masses , intensities,  ppm = ppm , z = z)
  }
  return(isotope_molecules_list)}

#' @title Filtering by validity
#' @description A list of annotation information of isotopes
#' @param isotope_molecules_list list of annotate information for isotopes
#' @usage isotope_valid(isotope_molecules_list)
isotope_valid <- function(isotope_molecules_list){
  isotope_valid_list <- NULL
  isotope_molecules_list[sapply(isotope_molecules_list, is.null)] <- NULL
  length_iso_v <- length(isotope_molecules_list)
  for (i in 1:length_iso_v) {
    if (isotope_molecules_list[[i]]$valid == "Valid") {
      isotope_valid_list[[i]] <- lapply(isotope_molecules_list[[i]], `[[`, 1)
    }
  }
  return(isotope_valid_list)
}


#' @title Creating dataframe from isotope_valid_list
#' @description Creating dataframe from isotope_valid_list
#' @param isotope_valid_list list of annotate information for isotopes
#' @usage valid_list_as_df(isotope_valid_list)
valid_list_as_df <- function(isotope_valid_list){
  isotope_valid_list[sapply(isotope_valid_list, is.null)] <- NULL
  rep <- length(isotope_valid_list)
  isotope_mass <- rep(NA, rep)
  exactmass <- rep(NA, rep)
  formula <- rep(NA, rep)
  score <- rep(NA, rep)
  charge <- rep(NA, rep)
  parity <- rep(NA, rep)
  DBE <- rep(NA, rep)
  valid <- rep(NA, rep)
  
  for (i in 1:length(isotope_valid_list)) {
    isotope_mass[i] <-  isotope_valid_list[[i]]$isotopes[1,2]
    exactmass[i] <-     isotope_valid_list[[i]]$exactmass
    formula[i] <-     isotope_valid_list[[i]]$formula
    score[i] <-     isotope_valid_list[[i]]$score
    charge[i] <-     isotope_valid_list[[i]]$charge
    parity[i] <-     isotope_valid_list[[i]]$parity
    DBE[i] <-     isotope_valid_list[[i]]$DBE
    valid[i] <-     isotope_valid_list[[i]]$valid
  }
  isotope_valid_df <- as.data.frame(cbind(valid,formula,exactmass,isotope_mass,score,charge,parity,DBE))
  return(isotope_valid_df)}



#' @title Mass spectrometer visualization
#' @description Gives a graphical representation of the mono-isotopic and isotopic peaks
#' @param final_df_vis dataframe of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  plotly plot_ly highlight_key add_markers add_segments  layout rangeslider ggplotly highlight
#' @usage mass_spec_plot(final_df_vis)
#' @importFrom dplyr %>%
mass_spec_plot <- function(final_df_vis){
  d <- highlight_key(final_df_vis, ~id )
  
  fig <- plot_ly(data = d, x = ~mz, y = ~intensity, color = ~Isotopic_status,colors = "Set1",
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
#' @param df Dataframe must contain columns c("mz", "intensity") with "mz" as first column
#' @param ppm Parts Per Million Tolerance
#' @param z charge
#' @param isotope_da difference in mass for monoisotopic molecule and isotope in Da
#' @usage isotope_tagging(df , ppm ,z , isotope_da)
#' @importFrom  data.table fintersect as.data.table
#' @importFrom dplyr mutate  filter select mutate case_when
#' @importFrom tibble tibble
#' @importFrom tidyselect everything
#' @importFrom sqldf sqldf
#' @export
isotope_tagging <- function(df, ppm, z, isotope_da){
  #-----------------------------------------------------------------------------------
  ## Start body part of function BAIT
  #keep original copy
  copy_df <- df 
  
  
  #creating new data frame with mean intensity column # note now complete case but will be replaced by imputation method
  df <- data.frame(mz = df[,1], intensity = rowMeans(df[,-1], na.rm = TRUE))
  
  #setting mz & intensity to vectors
  mz_vector <-  df$mz
  
  #creating new df with mass_error , lower bound and upperbound
  df <- mutate(df,
               mass_error = ppm_to_dalton(mz, ppm),
               min_error  = mz - mass_error,
               max_error  = mz + mass_error)
  #assigning minimum error and max error
  min_error_vec <- df$min_error
  max_error_vec <- df$max_error
  
  #-------------------------------------------------------------------------------------
  # Annotation
  
  # Annotation by mass difference 
  # computing maximum difference and minimum difference
  max_diff_mat <- abs(outer(max_error_vec, min_error_vec, `-`))
  min_diff_mat <- abs(outer(min_error_vec, max_error_vec, `-`))
  
  #filtering + indexing max and min difference
  index_max_mz <- which(max_diff_mat >= isotope_da, 
                        arr.ind = TRUE )
  index_min_mz <- which(min_diff_mat <= isotope_da & min_diff_mat > 0,
                        arr.ind = TRUE )
  
  
  index_dt_min <- as.data.table(index_min_mz )
  index_dt_max <- as.data.table(index_max_mz)
  
  index_dt_matches <- fintersect(index_dt_min, index_dt_max , all = TRUE) 
  
  #making new df of the isotopes 
  iso_df <- tibble(mol_ion = mz_vector[index_dt_matches[,col]],
                   isotope = mz_vector[index_dt_matches[,row]],) 
  
  # Adding intensities
  adducts_matches_int <- add_intensities(df, iso_df)
  
  # filter for isotope elements
  #adducts_isotope <- filter( adducts_matches_int,  adducts_matches_int$matches %in% c("13C"))
  
  # decompose isotopes
  isotope_molecules_list <- isotope_molecules(adducts_matches_int, ppm, z)
  
  # filtering keeping only valid results
  isotope_valid_list <- isotope_valid(isotope_molecules_list)
  
  # transforming isotope_valid_mist into dataframe
  isotope_valid_df <- valid_list_as_df(isotope_valid_list)
  
  # Filter for formula with "C" in it
  isotope_valid_df <- isotope_valid_df[grep("^C",isotope_valid_df$formula),]
  
  # Change character vectors to numeric
  isotope_valid_df$exactmass <- as.numeric(isotope_valid_df$exactmass)
  isotope_valid_df$isotope_mass <- as.numeric(isotope_valid_df$isotope_mass)
  
  # Merge together by condition between error values
  merge_df <- sqldf("select * from isotope_valid_df left join df
      on (isotope_valid_df.isotope_mass > df.min_error and isotope_valid_df.isotope_mass <= df.max_error) ")
  merge_df <- na.omit(merge_df)
  merge_dff <- sqldf("select * from isotope_valid_df inner join df on (isotope_valid_df.exactmass > df.min_error and isotope_valid_df.exactmass <= df.max_error) ")
  merge_dff <- na.omit(merge_dff)
  
  final <- merge(x = merge_df, y = merge_dff, by ="formula")
  final <- na.omit(final)
  final$id <- paste0('id-', 1:nrow(final))
  final$Isotopic_status <- paste0("Molecular_ion")
  final$iso_ratio <- final$intensity.x/final$intensity.y
  final_m <- final[,c(1,27,15,21,24,25,22,16,9,12,13,10,28,26,5)]
  colnames(final_m) <- c("fomula","Isotopic_status","exactmass","mz","min_error","max_error","intensity",
                         "isotope_mass","isotope_mz","iso_min_error","iso_max_error","intensity_iso","intensity_ratio","id","score")
  
  final$Isotopic_status <- paste0("c13_isotope")
  final_i <- final[,c(1,27,15,9,11,12,10,16,9,12,13,10,28,26,5)]
  colnames(final_i) <- c("fomula","Isotopic_status","exactmass","mz","min_error","max_error","intensity",
                         "isotope_mass","isotope_mz","iso_min_error","iso_max_error","intensity_iso","intensity_ratio","id","score")
  final_c <- rbind(final_m, final_i)
  final_c <- final_c[order(final_c$mz),]
  
  
  #Adding adducts to dataframe
  final_f <- merge(copy_df,final_c ,by="mz", all= TRUE)
  final_f <- final_f[,-c(1:ncol(copy_df))]
  #plotting
  #print(mass_spec_plot(final_c))
  
  rownames <- rownames(copy_df)
  #return data 
  return_df <- cbind(final_f, copy_df)
  rownames(return_df) <- rownames
  return(return_df)
}

