###############################################
#               sumR                          #
# Topic : Detection of isotopes and adducts   #
# Merys Valdez, February  2022,  LACDR        #
###############################################


# Content of script
# -----------------

# Part 1 : Functions to use in isotope_tagging 
#(1) ppm_to_dalton function
#(2) add_intensities function
#(3) isotope_molecules function 
#(4) isotope_valid function 
#(5) valid_list_as_df function 
#(6) isotope_finding function 
#(7) adduct_finding function 
#(8) big_merge function
#(9)  bar_plot   function  
#(10) mass_spec_plot function

# Part 2 : The isotope_tagging function
#(11) isotope_tagging function






#-------------------------------------------------------------------------------
##############################################
#               Part 1 :                     #
# Functions to use in isotope_tagging        #
##############################################

#-----------------------------
#(1) ppm_to_dalton function  - 
#-----------------------------
#' @title Dalton tolerance calculation
#' @description This function converts the Parts Per Million(PPM)tolerance into Dalton tolerance
#' @param mass Mass of interest in m/z
#' @param ppm  Parts Per Million tolerance
#' @usage ppm_to_dalton(mass, ppm)
ppm_to_dalton <- function(mass, ppm) {
  da_error  <- (mass * ppm) / 1e6
  return(da_error)}



#------------------------------
#(2) add_intensities function -
#------------------------------
#' @title Adding intensities
#' @description Adds an intensity column for mono-isotopic ion and isotope
#' @param df Data frame
#' @param iso_df Data frame with containing possible mono-isotopic ion and isotope based on solely mass differences
#' @usage add_intensities(df, iso_df)
#' @importFrom dplyr filter 
add_intensities <- function(df, iso_df){
  # Adding intensity for mono-isotopic ion
  add_mol <- merge(filter(df, df$mz %in% iso_df$mol_ion),iso_df, by.x = "mz",
                   by.y = "mol_ion", all.x = TRUE)
  add_mol <- add_mol[,c(1,2)]
  add_mol$intensity_mol <- add_mol$intensity
  add_mol$mz_mol <- add_mol$mz
  add_mol$intensity <- NULL
  add_mol$mz <- NULL
  
  # Adding intensity for isotope ion
  add_iso <- merge(filter(df, df$mz %in% iso_df$isotope),iso_df, by.x = "mz",
                   by.y = "isotope", all.x = TRUE)
  add_iso <- add_iso[,c(1,2)]
  add_iso$intensity_iso <- add_iso$intensity
  add_iso$mz_iso <- add_iso$mz
  add_iso$intensity <- NULL
  add_iso$mz <-NULL
  
  # Binding mono- and isotope masses and intensities
  iso_df_fin <- cbind(add_mol,add_iso)
  
  
  return(iso_df_fin)
}



#--------------------------------
#(3) isotope_molecules function -
#--------------------------------
#' @title List of isotopes with annotation information
#' @description Creates a list of isotopes with annotation information
#' @param iso_df_fin Data frame with containing mass and intensity of possible mono-isotopic ion and isotope pairs
#' @param ppm Parts Per Million tolerance
#' @param z Charge
#' @usage isotope_molecules(iso_df_fin, ppm , z)
#' @importFrom Rdisop decomposeIsotopes initializeElements
isotope_molecules <- function(iso_df_fin, ppm ,z){
  # Creating list for return
  isotope_molecules_list <- NULL
  
  # Initialize Elements
  element_vector <- c("C","H","N","O","P","S","Cl")
  element_vector  <- initializeElements(element_vector) 
  
  # Decompose isotopes
  length_iso <- nrow(iso_df_fin)
  for (i in 1:length_iso) {
    masses <- c(iso_df_fin$mz_mol[i], iso_df_fin$mz_iso[i])
    intensities <- c(100, ((iso_df_fin$intensity_iso[i] / iso_df_fin$intensity_mol[i]) * 100))
    isotope_molecules_list[[i]] <- decomposeIsotopes(masses, 
                                                     intensities, 
                                                     ppm = ppm, 
                                                     z = z, 
                                                     elements = element_vector)
  }
  return(isotope_molecules_list)}



#----------------------------
#(4) isotope_valid function -
#----------------------------
#' @title Filtering isotopes by validity
#' @description Filtering a list of possible isotopes based on validity
#' @param isotope_molecules_list list of possible isotopes containing annotate information 
#' @importFrom  purrr map
#' @usage isotope_valid(isotope_molecules_list)
isotope_valid <- function(isotope_molecules_list){
  # Creating list for return
  isotope_valid_list <- NULL
  
  # Removing missing values from list
  isotope_molecules_list[sapply(isotope_molecules_list, is.null)] <- NULL
  
  # Setting the lenght, to loop over
  length_iso_v <- length(isotope_molecules_list)
  
  # Loop searching/filtering for valid compounds
  for (i in 1:length_iso_v) {
    if (isotope_molecules_list[[i]]$valid == "Valid") {
      isotope_valid_list[[i]] <- map(isotope_molecules_list[[i]],1)
    }
  }
  return(isotope_valid_list)
}



#-------------------------------
#(5) valid_list_as_df function -
#-------------------------------
#' @title Converting isotope_valid_list  
#' @description Converting isotope_valid_list to a data frame
#' @param isotope_valid_list list of valid annotate information concerning the compounds
#' @usage valid_list_as_df(isotope_valid_list)
valid_list_as_df <- function(isotope_valid_list){
  
  # Removing missing values from list
  isotope_valid_list[sapply(isotope_valid_list, is.null)] <- NULL
  
  # Creating vectors/variables for each element in the list
  rep <- length(isotope_valid_list)
  isotope_mass <- rep(NA, rep)
  exactmass <- rep(NA, rep)
  formula <- rep(NA, rep)
  score <- rep(NA, rep)
  charge <- rep(NA, rep)
  parity <- rep(NA, rep)
  DBE <- rep(NA, rep)
  valid <- rep(NA, rep)
  
  # Looping for each list in the list to assign it to the matching vector
  for (i in 1:length(isotope_valid_list)) {
       isotope_mass[i] <-  isotope_valid_list[[i]]$isotopes[1,2]
       exactmass[i] <-     isotope_valid_list[[i]]$exactmass
       formula[i] <-       isotope_valid_list[[i]]$formula
       score[i] <-         isotope_valid_list[[i]]$score
       charge[i] <-        isotope_valid_list[[i]]$charge
       parity[i] <-        isotope_valid_list[[i]]$parity
       DBE[i] <-           isotope_valid_list[[i]]$DBE
       valid[i] <-         isotope_valid_list[[i]]$valid
  }
  
  # Creating data frame by combining the vectors annotation variables
  isotope_valid_df <- as.data.frame(cbind(valid, 
                                          formula,
                                          exactmass, 
                                          isotope_mass,
                                          score, 
                                          charge,
                                          parity,
                                          DBE))
  return(isotope_valid_df)}



#-------------------------------
#(6) isotope_finding function  -
#-------------------------------
#' @title Isotope finding
#' @description Searches for isotopes 
#' @param df Input data frame, with min_error and max_error columns added
#' @param mz_vector vector of mz values
#' @param isotope_df Data frame containing isotopic elements of interest, PSE notation and mass difference compared to mono-isotope
#' @param Elements Vector with selected isotopic elements of interest 
#' @param ppm Parts Per Million tolerance
#' @param z Charge
#' @param max_diff_mat Matrix with maximum difference between mz pairs
#' @param min_diff_mat Matrix with maximum difference between mz pairs
#' @importFrom tibble tibble
#' @importFrom data.table as.data.table fintersect
#' @importFrom dplyr filter
#' @usage isotope_finding(df, mz_vector,isotope_df, Elements, ppm, z , max_diff_mat,min_diff_mat)
isotope_finding <- function(df, mz_vector, isotope_df, Elements, ppm, z, max_diff_mat, min_diff_mat){
  
  # Filtering for isotopic elements of interest
  isotope_filter <- filter(isotope_df, isotope_df$Element %in% Elements)
  
  # Creating data frame for return
  results_valid <- NULL
  
  # Looping for each isotopic elements of interest
  for (i in 1:length(isotope_filter$isotope_da)){
    
    # Searching for mass differences that fall within tolerance interval for each possible mono- and isotopic pair
    index_max_mz <- which(max_diff_mat >= isotope_filter$isotope_da[i], 
                          arr.ind = TRUE )
    index_min_mz <- which(min_diff_mat <= isotope_filter$isotope_da[i]& 
                          min_diff_mat > 0,
                          arr.ind = TRUE )
    
    # Transforming to data frame to data table to improve speed of merging
    index_dt_min <- as.data.table(index_min_mz)
    index_dt_max <- as.data.table(index_max_mz)
    
    # Selecting mass differences that fall within tolerance interval for each possible mono- and isotopic pair
    index_dt_matches <- fintersect(index_dt_min, index_dt_max, all = TRUE) 
    
    # Creating a data frame for mono- and isotopic ions 
    iso_df <- tibble(mol_ion = mz_vector[index_dt_matches[,col]],
                     isotope = mz_vector[index_dt_matches[,row]],) 
    
    # Adding intensity values for each mz value
    adducts_matches_int <- add_intensities(df, iso_df)
    
    # Decompose isotopes
    isotope_molecules_list <- isotope_molecules(adducts_matches_int, ppm, z)
    
    # Filtering keeping only valid compounds
    isotope_valid_list <- isotope_valid(isotope_molecules_list)
    
    # Transforming isotope_valid_list into a data frame
    isotope_valid_df <- valid_list_as_df(isotope_valid_list)
    
    # Checking if formula contains "X"- element of interest
    isotope_valid_df <- isotope_valid_df[grep(isotope_df$Element_notation[i],
                                              isotope_valid_df$formula),]
    
    # Transforming character vectors to numeric
    isotope_valid_df$exactmass <- as.numeric(isotope_valid_df$exactmass)
    isotope_valid_df$isotope_mass <- as.numeric(isotope_valid_df$isotope_mass)
    
    # Iterative row-binding for each "X"- element of interest
    results_valid <- rbind(results_valid,isotope_valid_df)
  }
  return(results_valid)
}



#-------------------------------
#(7) adduct_finding  function  -
#-------------------------------
#' @title Adduct finding
#' @description Searches for adducts based on mass difference
#' @param adduct_df Data frame containing adduct information, name and mass columns
#' @param max_diff_mat Matrix with maximum difference between mz pairs
#' @param min_diff_mat Matrix with minumum difference between mz pairs
#' @param mz_vector vector of mz values
#' @importFrom  data.table as.data.table fintersect 
#' @importFrom  tibble tibble
#' @usage adduct_finding(adduct_df, max_diff_mat, min_diff_mat, mz_vector)
adduct_finding <- function(adduct_df, max_diff_mat, min_diff_mat, mz_vector){
  # Creating data frame for return
  results_valid <- NULL
  
  # Looping to search for each mass of an adduct
  for (i in 1:length(adduct_df$mass)){
       index_max_mz <- which(max_diff_mat >= adduct_df$mass[i], 
                          arr.ind = TRUE )
       index_min_mz <- which(min_diff_mat <= adduct_df$mass[i] & min_diff_mat > 0,
                          arr.ind = TRUE )
    
       # Transforming from data frame to data table to optimize speed of merging 
       index_dt_min <- as.data.table(index_min_mz)
       index_dt_max <- as.data.table(index_max_mz)
    
       index_dt_matches <- fintersect(index_dt_min, index_dt_max , all = TRUE) 
    
       # Creating new dataframe for adducts 
       iso_df <- tibble(mol_ion = mz_vector[index_dt_matches[,col]],
                     isotope = mz_vector[index_dt_matches[,row]],
                     adduct  = adduct_df$name[i]) 
       results_valid <- rbind(results_valid, iso_df)
  }
  return(results_valid)
}



#-------------------------------
#(8) big_merge function        -
#-------------------------------
#' @title Big Merging function
#' @description Merging all different isotope data frames into one
#' @param isotope_valid_df Data frame of valid annotate information for isotopes
#' @param df Input data frame, with min_error and max_error columns added
#' @param copy_df A copy of the input data frame
#' @importFrom  sqldf sqldf
#' @importFrom dplyr %>% 
#' @importFrom tibble add_column
#' @usage big_merge(isotope_valid_df, df, copy_df)
big_merge <- function(isotope_valid_df, df, copy_df){
  
  # Merge based upon the condition of falling within min and max error mass for isotopes
  merge_df <- sqldf("select * from isotope_valid_df left join df on 
                    (isotope_valid_df.isotope_mass > df.min_error and isotope_valid_df.isotope_mass <= df.max_error) ")
  # Merge based upon the condition of falling within min and max error mass for isotopes
  merge_dff <- sqldf("select * from isotope_valid_df inner join df on 
                     (isotope_valid_df.exactmass > df.min_error and isotope_valid_df.exactmass <= df.max_error) ")
  
  # Create final data frame containing annotation information for mono-isotope
  final <- merge(x = merge_df, y = merge_dff, by = "formula")
  final <- na.omit(final)
  final$id <- paste0('id-', 1:nrow(final))
  final$isotopic_status <- paste0("Molecular_ion")
  final$iso_ratio <- round((final$intensity.x/final$intensity.y),2)
  final_m <- final[,c(1, 15, 27, 21, 24, 25, 22, 16, 9, 12, 13, 10, 28, 26, 5)]
  colnames(final_m) <- c("fomula", "exactmass_compound", "isotopic_status", 
                         "mz", "min_error", "max_error", "intensity",
                         "isotope_mass", "isotope_mz", "iso_min_error", 
                         "iso_max_error", "intensity_iso", "intensity_ratio",
                         "id", "score")
  
  # Create final data frame containing annotation information for isotope
  
  # Assigning status depending if "X"-element found in formula
  final$isotopic_status  <- ifelse(grepl("Cl", final$formula) == TRUE ,"Cl37-isotope", "C13-isotope" )
  
  final_i <- final[,c(1, 15, 27, 9, 11, 12, 10, 16, 9, 12, 13, 10, 28, 26, 5)]
  colnames(final_i) <- c("fomula", "exactmass_compound", "isotopic_status", 
                         "mz", "min_error", "max_error", "intensity",
                         "isotope_mass", "isotope_mz", "iso_min_error", 
                         "iso_max_error", "intensity_iso", "intensity_ratio",
                         "id", "score")
  final_c <- rbind(final_m, final_i)
  final_c <- final_c[order(final_c$mz),]
  
  
  # Correcting for duplicates
  final_f <- merge(copy_df, final_c ,by = "mz", all = TRUE)
  final_f <- final_f %>% distinct(mz, .keep_all = TRUE)
  final_f <- final_f[,-c(1:ncol(copy_df))]
  
  # Setting mz in line
  result <- add_column(final_f, mz =copy_df$mz, .after = 3)
  
  return(result)
}

#-------------------------------
#(9)       bar_plot   function  -
#-------------------------------
#' @title Bar plot visualization
#' @description Displays a bar plot for the different isotopic status 
#' @param final_df data frame of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  ggplot2 ggplot geom_bar scale_fill_brewer xlab geom_text ylab
#' @usage bar_plot(final_df)
#' @importFrom dplyr %>% rename
bar_plot <- function(final_df){
  
# Filtering for only available isotopic status 
vis_df <- final_df[!is.na(final_df$isotopic_status), ]
vis_df$isotopic_status <- factor(vis_df$isotopic_status)

# Creating data frame for barplot
count_df <- as.data.frame(table(vis_df$isotopic_status))   
count_df  <- count_df %>% 
  rename(
    isotopic_status = Var1,
    count = Freq
  )
# Bar plot
p <- ggplot(count_df, aes(x=isotopic_status,y=count, fill=isotopic_status))+
  geom_bar(stat="identity", width=0.7)+
  scale_fill_brewer(palette="Dark2") +
  xlab(label = "Isotopic status") +
  geom_text(aes(label = count), vjust = -1) +
  ylab(label="Count") 
return(p)
}
#-------------------------------
#(10) mass_spec_plot  function  -
#-------------------------------
#' @title Mass spectrometry visualization
#' @description Gives a graphical representation of the mono-isotopic and isotopic peaks
#' @param final_df dataframe of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  plotly plot_ly highlight_key add_markers add_segments  layout rangeslider ggplotly highlight
#' @usage mass_spec_plot(final_df)
#' @importFrom dplyr %>%
mass_spec_plot <- function(final_df){
  vis_df <- final_df[!is.na(final_df$isotopic_status), ]
  d <- highlight_key(vis_df , ~id)
  
  fig <- plot_ly(data = d, 
                 x = ~mz, 
                 y = ~intensity, 
                 color = ~isotopic_status,
                 colors = "Dark2",
                 hoverinfo = 'text',
                 text = ~paste("<b>m/z</b>: ",round(mz,3), '<br><b>id</b>:', id, '<br><b>intensity ratio</b>:', intensity_ratio))
  fig <- fig %>%
    add_markers() %>%
    add_segments(
      data = vis_df , x = ~mz, xend = ~mz,
      y = 0, yend = ~intensity, color = I("grey"), showlegend = FALSE)  %>% 
    layout(
        title = "Mass spectrometry of detected isotopes")
      
  
  gg <- ggplotly(fig)
  gg <- highlight(gg, dynamic = TRUE, off ='plotly_doubleclick' )
  gg  
  
  
  return(gg)
}



#-------------------------------------------------------------------------------

##############################################
#               Part 2 :                     #
#    The isotope_tagging function            #
##############################################


#---------------------------------
#(11) isotope_tagging  function  -
#---------------------------------
#' @title Isotope tagging
#' @description Tags isotopes based on intensity ratio and difference in mass between mono- and isotopic ion
#' @param df Input data frame which must contain columns c("mz", "intensity") with "mz" as first column
#' @param ppm Parts Per Million tolerance
#' @param z Charge
#' @usage isotope_tagging(df, ppm ,z)
#' @importFrom data.table fintersect as.data.table
#' @importFrom dplyr %>% mutate  distinct
#' @export
isotope_tagging <- function(df, ppm, z){
  #-----------------------------------------------------------------------------
  # Keep original copy
  copy_df <- df 
  
  
  #-----------------------------------------------------------------------------
  # Creating isotope and adduct information data frames
  
  # Isotope data frame information
  Element_notation <- c("C","Cl")
  Element <- c("C13","Cl37")
  isotope_da <- c(1.0034 , 1.9970)
  isotope_df <- as.data.frame(cbind(Element, Element_notation, isotope_da))
  isotope_df$isotope_da <- as.numeric(isotope_df$isotope_da)
  
  # Adduct data frame information
  name <- c("Na adduct"," K adduct ")
  mass <- c(21.98194 , 37.95589)
  adduct_df <- as.data.frame(cbind(name, mass))
  adduct_df$mass <- as.numeric(adduct_df$mass)
  
  
  # Start body of isotope_tagging function
  #-----------------------------------------------------------------------------
  # Creating new data frame with mean intensity column # note now complete case but will be replaced by imputation method
  df <- data.frame(mz = df[,1], intensity = rowMeans(df[,-1], na.rm = TRUE))
  
  # Setting mz variables to vectors
  mz_vector <-  df$mz
  
  # Creating new data frame with mass_error , lower bound and upper bound
  df <- mutate(df,
               mass_error = ppm_to_dalton(mz, ppm),
               min_error  = mz - mass_error,
               max_error  = mz + mass_error)
  
  # Assigning minimum error and max error
  min_error_vec <- df$min_error
  max_error_vec <- df$max_error
  
  
  # Computing maximum difference and minimum difference between mz values
  max_diff_mat <- abs(outer(max_error_vec, min_error_vec, `-`))
  min_diff_mat <- abs(outer(min_error_vec, max_error_vec, `-`))
  
  #-----------------------------------------------------------------------------
  # Searching for valid isotopes of interest
  isotope_valid_df <- isotope_finding(df, mz_vector,isotope_df,Elements=c("C13","Cl37"), ppm, z , max_diff_mat,min_diff_mat)
  
  
  
  # Big merge and printing the number of found isotopes
  final_f  <- big_merge(isotope_valid_df, df, copy_df)
  print(paste("Number of isotopes found:",length(unique(final_f$id))))
  
  
  #-----------------------------------------------------------------------------
  # Searching for adducts
  adduct_found <- adduct_finding(adduct_df, max_diff_mat, min_diff_mat, mz_vector)
  test_adduct <- merge(final_f, adduct_found[,2:3], by.x ="mz" , by.y = "isotope", all.x = TRUE)
  final_df <- test_adduct %>% distinct(mz, .keep_all = TRUE)
  
  
  
  #-----------------------------------------------------------------------------
  # Returning original features with data 
  
  # Extracting feature names
  rownames <- rownames(copy_df)
  
  # Return data frame
  return_df <- cbind(final_df, copy_df)
  rownames(return_df) <- rownames
  
  
  
  
  #-----------------------------------------------------------------------------
  # Plotting
   print(bar_plot(final_df))
   print(mass_spec_plot(final_df))
  
  

  
  return(return_df)
}



  
