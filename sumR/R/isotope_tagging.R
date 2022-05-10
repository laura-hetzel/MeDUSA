###############################################
#               sumR                          #
# Topic : Detection of isotopes and adducts   #
# Merys Valdez, February  2022,  LACDR        #
###############################################


# Content of script
# -----------------

# Part 1 : Functions to use in isotope_tagging
# (1) ppmToDalton function
# (2) addIntensities function
# (3) isotopeMolecules function
# (4) isotopeValid function
# (5) validListAsDf function
# (6) isotopeFinding function
# (7) adductFinding function
# (8) bigMerge function
# (9) barPlot   function
# (10)massSpecPlot function

# Part 2 : The isotopeTagging function
# (11) isotopeTagging function






#-------------------------------------------------------------------------------
##############################################
#               Part 1 :                     #
# Functions to use in isotope_tagging        #
##############################################

#-----------------------------
# (1) ppmToDalton function  -
#-----------------------------
#' @title Dalton tolerance calculation
#' @description This function converts the parts per million(ppm)tolerance into Dalton tolerance unit
#' @param mass A vector containing mass of interest in m/z
#' @param ppm  An integer, defining parts per million (ppm) for tolerance (default = 5)
#' @usage ppmToDalton(mass, ppm)
ppmToDalton <- function(mass, ppm = 5) {
  da_error <- (mass * ppm) / 1e6
  return(da_error)
}



#------------------------------
# (2) addIntensities function -
#------------------------------
#' @title Adding intensities
#' @description Adds an intensity column for mono-isotopic ion and isotope ions
#' @param data data.frame, containing mz and intensity information
#' @param iso_df data.frame that contains possible mono-isotopic ion and isotope based on solely mass differences
#' @importFrom dplyr filter
#' @usage addIntensities(data, iso_df)
addIntensities <- function(data, iso_df) {
  # Adding intensity for mono-isotopic ion
  add_mol <- merge(dplyr::filter(data, data$mz %in% iso_df$mol_ion), iso_df,
    by.x = "mz",
    by.y = "mol_ion", all.x = TRUE
  )
  add_mol <- add_mol[, c(1, 2)]
  add_mol <- data.frame(intensity_mol = add_mol$intensity, mz_mol = add_mol$mz)


  # Adding intensity for isotope ion
  add_iso <- merge(dplyr::filter(data, data$mz %in% iso_df$isotope), iso_df,
    by.x = "mz",
    by.y = "isotope", all.x = TRUE
  )
  add_iso <- add_iso[, c(1, 2)]


  # Binding mono- and isotope masses and intensities into data frame
  iso_df_fin <- data.frame(
    intensity_mol = add_mol$intensity, mz_mol = add_mol$mz,
    intensity_iso = add_iso$intensity, mz_iso = add_iso$mz
  )


  return(iso_df_fin)
}



#--------------------------------
# (3) isotopeMolecules function -
#--------------------------------
#' @title List of isotopes with annotation information
#' @description Creates a list of isotopes with annotation information
#' @param adducts_matches_int data.frame with containing mass and intensity of possible mono-isotopic ion and isotope pairs
#' @param ppm An integer, defining parts per million (ppm) for tolerance (default = 5)
#' @param z An integer, defining charge z of m/z peaks for calculation of real mass. 0 is for auto-detection (default = 0)
#' @param Elements vector, with selected isotopic elements of interest (default = c("C13")), options : c("C13","Cl37")
#' @importFrom Rdisop decomposeIsotopes initializeElements
#' @usage isotopeMolecules(adducts_matches_int, Elements = c("C13"), ppm = 5 , z = 0)
isotopeMolecules <- function(adducts_matches_int, Elements = c("C13"), ppm = 5, z = 0) {

  # Check isotopic selection
  if (length(Elements) == 2) {
    element_vector <- c("C", "H", "N", "O", "P", "S", "Cl")
  } else {
    element_vector <- c("C", "H", "N", "O", "P", "S")
  }

  # Initialize Elements
  element_vector <- initializeElements(element_vector)

  # Decompose isotopes
  isotope_molecules_list <- lapply(1:nrow(adducts_matches_int), function(i) {
    masses <- c(adducts_matches_int$mz_mol[i], adducts_matches_int$mz_iso[i])
    intensities <- c(100, ((adducts_matches_int$intensity_iso[i] / adducts_matches_int$intensity_mol[i]) * 100))
    decomposeIsotopes(masses, intensities, ppm = ppm, z = z, elements = element_vector)
  })

  return(isotope_molecules_list)
}



#----------------------------
# (4) isotopeValid function -
#----------------------------
#' @title Filtering isotopes by validity
#' @description Filtering a list of possible isotopes based on validity
#' @param isotope_molecules_list A list, of possible isotopes containing annotate information
#' @usage isotopeValid(isotope_molecules_list)
isotopeValid <- function(isotope_molecules_list) {
  # Removing missing values from list
  isotope_molecules_list[sapply(isotope_molecules_list, is.null)] <- NULL

  # Loop searching/filtering for valid compounds
  isotope_valid_list <- lapply(1:length(isotope_molecules_list), function(i) {
    if (isotope_molecules_list[[i]]$valid[1] == "Valid") {
      lapply(isotope_molecules_list[[i]], `[[`, 1)
    }
  })
  return(isotope_valid_list)
}


#-------------------------------
# (5) validListAsDf function -
#-------------------------------
#' @title Converting isotope_valid_list
#' @description Converting isotope_valid_list to a data frame
#' @param isotope_valid_list A list, of valid annotate information concerning the compounds
#' @usage validListAsDf(isotope_valid_list)
validListAsDf <- function(isotope_valid_list) {

  # Removing missing values from list
  isotope_valid_list[sapply(isotope_valid_list, is.null)] <- NULL

  # Transforming list into data frame
  isotope_valid_df <- do.call(rbind, lapply(isotope_valid_list, function(i) {
    data.frame(
      valid = i$valid,
      formula = i$formula,
      exactmass = i$exactmass,
      isotope_mass = i$isotopes[1, 2],
      score = i$score,
      charge = i$charge,
      parity = i$parity,
      DBE = i$DBE
    )
  }))

  return(isotope_valid_df)
}


#-------------------------------
# (6) isotopeFinding function  -
#-------------------------------
#' @title Isotope finding
#' @description Searches for isotopes
#' @param data data.frame, with min_error and max_error columns added
#' @param mz_vector vector of mz values
#' @param isotope_df data.frame containing isotopic elements of interest, PSE notation and mass difference compared to mono-isotope
#' @param Elements vector, with selected isotopic elements of interest (default = c("C13")), options : c("C13","Cl37")
#' @param ppm An integer, defining parts per million (ppm) for tolerance (default = 5)
#' @param z An integer, defining charge z of m/z peaks for calculation of real mass. 0 is for auto-detection (default = 0)
#' @param max_diff_mat matrix, with maximum difference between mz pairs
#' @param min_diff_mat matrix, with maximum difference between mz pairs
#' @importFrom data.table fintersect as.data.table
#' @importFrom dplyr filter
#' @usage isotopeFinding(data, mz_vector, isotope_df, Elements = c("C13"), ppm = 5, z = 0, max_diff_mat, min_diff_mat)
isotopeFinding <- function(data, mz_vector, isotope_df, Elements = c("C13"), ppm = 5, z = 0, max_diff_mat, min_diff_mat) {
  # Setting Elements of interest

  # Filtering for isotopic elements of interest
  isotope_filter <- dplyr::filter(isotope_df, isotope_df$Element %in% Elements)

  # Looping for each isotopic elements of interest
  results_valid <- do.call(rbind, lapply(1:length(isotope_filter$isotope_da), function(i) {

    # Searching for mass differences that fall within tolerance interval for each possible mono- and isotopic pair
    index_max_mz <- which(max_diff_mat >= isotope_filter$isotope_da[i],
      arr.ind = TRUE
    )
    index_min_mz <- which(min_diff_mat <= isotope_filter$isotope_da[i] &
      min_diff_mat > 0,
    arr.ind = TRUE
    )

    # Transforming to data frame to data table to improve speed of merging
    index_dt_min <- as.data.table(index_min_mz)
    index_dt_max <- as.data.table(index_max_mz)

    # Selecting mass differences that fall within tolerance interval for each possible mono- and isotopic pair
    index_dt_matches <- fintersect(index_dt_min, index_dt_max, all = TRUE)

    # Creating a data frame for mono- and isotopic ions
    iso_df <- data.frame(
      mol_ion = mz_vector[index_dt_matches[, col]],
      isotope = mz_vector[index_dt_matches[, row]]
    )

    # Adding intensity values for each mz value
    adducts_matches_int <- addIntensities(data, iso_df)

    # Decompose isotopes
    isotope_molecules_list <- isotopeMolecules(adducts_matches_int, Elements, ppm, z)

    # Filtering keeping only valid compounds
    isotope_valid_list <- isotopeValid(isotope_molecules_list)

    # Transforming isotope_valid_list into a data frame
    isotope_valid_df <- validListAsDf(isotope_valid_list)

    # Checking if formula contains "X"- element of interest
    isotope_valid_df <- isotope_valid_df[grep(isotope_df$Element_notation[i], isotope_valid_df$formula), ]

    # Iterative row-binding for each "X"- element of interest
    return(isotope_valid_df)
  }))
  return(results_valid)
}



#-------------------------------
# (7) adductFinding  function  -
#-------------------------------
#' @title Adduct finding
#' @description Searches for adducts based on mass differences
#' @param adduct_df data.frame containing adduct information, name and mass columns
#' @param max_diff_mat matrix with maximum difference between mz pairs
#' @param min_diff_mat matrix with minumum difference between mz pairs
#' @param mz_vector vector of mz values
#' @importFrom  data.table as.data.table fintersect
#' @usage adductFinding(adduct_df, max_diff_mat, min_diff_mat, mz_vector)
adductFinding <- function(adduct_df, max_diff_mat, min_diff_mat, mz_vector) {

  # Looping to search for each mass of an adduct and row binding them
  results_valid <- do.call(rbind, lapply(1:length(adduct_df$mass), function(i) {
    index_max_mz <- which(max_diff_mat >= adduct_df$mass[i],
      arr.ind = TRUE
    )
    index_min_mz <- which(min_diff_mat <= adduct_df$mass[i] & min_diff_mat > 0,
      arr.ind = TRUE
    )

    # Transforming from data frame to data table to optimize speed of merging
    index_dt_min <- as.data.table(index_min_mz)
    index_dt_max <- as.data.table(index_max_mz)

    # merging found matches
    index_dt_matches <- fintersect(index_dt_min, index_dt_max, all = TRUE)

    # Creating new data frame for adducts
    adduct_df <- data.frame(
      compound = mz_vector[index_dt_matches[, col]],
      adduct = mz_vector[index_dt_matches[, row]],
      adducts = adduct_df$name[i]
    )
    return(adduct_df)
  }))
  return(results_valid)
}



#-------------------------------
# (8) bigMerge function        -
#-------------------------------
#' @title Big Merging function
#' @description Merging all different isotope data frames into one
#' @param isotope_valid_df a data.frame of valid annotate information for isotopes
#' @param data the input data.frame, with min_error and max_error columns bind
#' @param copy_df a data.frame copy of the input data frame
#' @importFrom  data.table setDT
#' @importFrom dplyr %>%
#' @usage bigMerge(isotope_valid_df, data, copy_df)
bigMerge <- function(isotope_valid_df, data, copy_df) {
  # merging isotope_valid_df witf data based on the condition that
  # the mass off compound or it's isotopes is between the min error and max error of mz

  # Setting data.frame to data.table
  setDT(isotope_valid_df)
  setDT(data)

  # merging for isotopes
  merge_isotope <- isotope_valid_df[data,
    on = .(isotope_mass >= min_error, isotope_mass < max_error), nomatch = 0,
    .(formula, exactmass, isotope_mass, score, charge, mz, min_error, max_error, intensity)
  ]
  # merging for mono - isotopes
  merge_mono <- isotope_valid_df[data,
    on = .(exactmass >= min_error, exactmass < max_error), nomatch = 0,
    .(formula, exactmass, isotope_mass, score, charge, mz, min_error, max_error, intensity)
  ]



  #-----------------------------------------------------------------------------
  # Create final data frame containing annotation information for mono-isotope
  merge_mono_iso <- merge(x = merge_isotope, y = merge_mono, by = "formula")

  # Adding id, isotopic status and intensity ratio
  merge_mono_iso$id <- paste0(1:nrow(merge_mono_iso))
  merge_mono_iso$isotopic_status <- paste0("Molecular_ion")
  merge_mono_iso$int_ratio <- round((merge_mono_iso$intensity.x / merge_mono_iso$intensity.y), 2)
  merge_mono_iso <- as.data.frame(merge_mono_iso)
  # Creating column names which need to be selected
  colnames_mono <- c(
    "formula", "exactmass.y", "isotopic_status", "mz.y", "min_error.y", "max_error.y", "intensity.y",
    "isotope_mass.y", "mz.x", "min_error.x", "max_error.x", "intensity.x", "int_ratio", "id", "score.y"
  )

  # Selecting columns for mono isotopic information
  final_mono <- merge_mono_iso[, colnames_mono]

  # Creating final column names
  colnames_mono_iso <- c(
    "formula", "exactmass_compound", "isotopic_status", "mz", "min_error", "max_error", "intensity",
    "isotope_mass", "isotope_mz", "iso_min_error", "iso_max_error", "intensity_iso", "intensity_ratio",
    "id", "score"
  )

  # Setting column names for mono isotopic information
  colnames(final_mono) <- colnames_mono_iso



  #-----------------------------------------------------------------------------
  # Create final data frame containing annotation information for isotope

  # Assigning status depending if "X"-element found in formula
  merge_mono_iso$isotopic_status <- ifelse(grepl("Cl", merge_mono_iso$formula) == TRUE, "Cl37-isotope", "C13-isotope")

  # Creating column names which need to be selected
  colnames_iso <- c(
    "formula", "exactmass.y", "isotopic_status", "mz.x", "min_error.x", "max_error.x", "intensity.x",
    "isotope_mass.x", "mz.x", "min_error.x", "max_error.x", "intensity.x", "int_ratio", "id", "score.y"
  )

  # Selecting columns for isotopic information
  final_iso <- merge_mono_iso[, colnames_iso]

  # Setting column names for isotopic information
  colnames(final_iso) <- colnames_mono_iso



  #-----------------------------------------------------------------------------
  # Row binding mono- and isotopic information
  final_mono_iso <- rbind(final_mono, final_iso)
  final_mono_iso <- final_mono_iso[order(final_mono_iso$mz), ]

  # Correcting for duplicates
  final_mono_iso <- merge(copy_df, final_mono_iso, by = "mz", all = TRUE)
  final_mono_iso <- final_mono_iso %>% distinct(mz, .keep_all = TRUE)
  final_mono_iso <- final_mono_iso[, -c(1:ncol(copy_df))]

  # Setting mz in line
  result <- cbind(final_mono_iso[, 1:3], mz = copy_df$mz, final_mono_iso[, -c(1:3)])
  result <- result %>% distinct(mz, .keep_all = TRUE)
  return(result)
}



#-------------------------------
# (9)      barPlot   function  -
#-------------------------------
#' @title Bar plot visualization
#' @description Displays a bar plot for the different types of isotopic status'
#' @param final_df a data.frame of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  ggplot2 ggplot geom_bar scale_fill_brewer xlab geom_text ylab
#' @importFrom dplyr rename %>%
#' @usage barPlot(final_df)
barPlot <- function(final_df) {

  # Filtering for only available isotopic status
  vis_df <- final_df[!is.na(final_df$isotopic_status), ]
  vis_df$isotopic_status <- factor(vis_df$isotopic_status)

  # Creating data frame for bar plot
  count_df <- as.data.frame(table(vis_df$isotopic_status))
  count_df <- count_df %>%
    rename(
      isotopic_status = Var1,
      count = Freq
    )
  # Bar plot
  bar_plot <- ggplot(count_df, aes(x = isotopic_status, y = count, fill = isotopic_status)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_brewer(palette = "Dark2", name = "Types of isotopic elements") +
    xlab(label = "Isotopic status") +
    ylab(label = "Count") +
    geom_text(aes(label = count), vjust = -1)
  return(bar_plot)
}



#-------------------------------
# (10) massSpecPlot  function  -
#-------------------------------
#' @title Mass spectrometry visualization
#' @description Gives a interactive graphical representation of the mono-isotopic and isotopic peaks
#' @param final_df a data.frame of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  plotly plot_ly highlight_key add_markers add_segments  layout rangeslider ggplotly highlight
#' @usage massSpecPlot(final_df)
#' @importFrom dplyr %>%
massSpecPlot <- function(final_df) {
  vis_df <- final_df[!is.na(final_df$isotopic_status), ]
  d <- highlight_key(vis_df, ~id)

  fig <- plot_ly(
    data = d,
    x = ~mz,
    y = ~intensity,
    color = ~isotopic_status,
    colors = "Dark2",
    hoverinfo = "text",
    text = ~ paste("<b>m/z</b>: ", round(mz, 3), "<br><b>id</b>:", id, "<br><b>intensity ratio</b>:", intensity_ratio)
  )
  fig <- fig %>%
    add_markers() %>%
    add_segments(data = vis_df, x = ~mz, xend = ~mz, y = 0, yend = ~intensity, color = I("grey"), showlegend = FALSE) %>%
    layout(title = "Mass spectrometry of detected isotopes")

  gg <- ggplotly(fig)
  gg <- highlight(gg, dynamic = TRUE, off = "plotly_doubleclick")
  gg


  return(gg)
}




#-------------------------------------------------------------------------------

##############################################
#               Part 2 :                     #
#    The isotope_tagging function            #
##############################################


#---------------------------------
# (11) isotopeTagging  function  -
#---------------------------------
#' @title Isotope tagging
#' @description Tags isotopes based on intensity ratio and difference in mass between mono- and isotopic ions
#' @param data data.frame,  which contain column named "mz", all other columns will be coerced as intensity
#' @param ppm An integer, defining parts per million (ppm) for tolerance (default = 5)
#' @param z An integer, defining charge z of m/z peaks for calculation of real mass. 0 is for auto-detection (default = 0)
#' @param Elements A vector containing the isotopic element of interest (default = c("C13"))
#' @param plot Logical, returns box plot for isotope and interactive mass spectrometry of detected isotopes (default = TRUE)
#' @importFrom dplyr distinct select %>%
#' @usage isotopeTagging(data, ppm = 5, Elements = c("C13"), z = 0, plot = TRUE)
#' @export
isotopeTagging <- function(exp, assay = "Area", ppm = 5, Elements = c("C13"), z = 0, plot = TRUE, filter = TRUE) {
  data <- as.data.frame(cbind(mz = rowData(exp)$mz, assay(exp, assay)))
  #-----------------------------------------------------------------------------
  # Check for input arguments
  if (is.null(data)) {
    stop("No data detected")
  }

  Check <- c("C13", "Cl37")
  if (!Elements %in% Check) {
    stop("Invalid 'Elements' value")
  }

  #-----------------------------------------------------------------------------
  # Keep original copy
  copy_df <- data



  #-----------------------------------------------------------------------------
  # Creating isotope and adduct information data frames

  # Isotope data frame information
  isotope_df <- data.frame(Element = c("C13", "Cl37"), Element_notation = c("C", "Cl"), isotope_da = c(1.0034, 1.9970))

  # Adduct data frame information
  adduct_df <- data.frame(name = c("Na adduct", "K adduct", "NH4 adduct"), mass = c(21.98194, 37.95589, 17.02660))



  # Start body of isotope_tagging function
  #-----------------------------------------------------------------------------
  # Creating new data frame with mean intensity column # note now complete case but will be replaced by imputation method
  df <- data.frame(mz = data[, c("mz")], intensity = rowMeans(select(data, -c("mz")), na.rm = TRUE))

  # Setting mz variables to vectors
  mz_vector <- as.numeric(df$mz)

  # Creating new data frame with mass_error , lower bound and upper bound
  df <- data.frame(mz = df$mz, intensity = df$intensity, mass_error = ppmToDalton(df$mz, ppm), min_error = df$mz - ppmToDalton(df$mz, ppm), max_error = df$mz + ppmToDalton(df$mz, ppm))

  # Assigning minimum error and max error
  min_error_vec <- df$min_error
  max_error_vec <- df$max_error

  # Computing maximum difference and minimum difference between mz values
  max_diff_mat <- abs(outer(max_error_vec, min_error_vec, `-`))
  min_diff_mat <- abs(outer(min_error_vec, max_error_vec, `-`))



  #-----------------------------------------------------------------------------
  # Searching for valid isotopes of interest
  isotope_valid_df <- isotopeFinding(df, mz_vector, isotope_df, Elements = Elements, ppm, z, max_diff_mat, min_diff_mat)


  # Big merge and printing the number of found isotopes
  big_merge_df <- bigMerge(isotope_valid_df, df, copy_df)
  print(paste("Number of isotopes found:", length(unique(big_merge_df$id))))



  #-----------------------------------------------------------------------------
  # Searching for adducts
  adduct_found <- adductFinding(adduct_df, max_diff_mat, min_diff_mat, mz_vector)
  test_adduct <- merge(big_merge_df, adduct_found[, 2:3], by.x = "mz", by.y = "adduct", all.x = TRUE)
  final_df <- test_adduct %>% distinct(mz, .keep_all = TRUE)


  #-----------------------------------------------------------------------------
  # Returning original features with data

  #-----------------------------------------------------------------------------
  # Plotting
  if (plot == TRUE) {
    print(barPlot(final_df))
    print(massSpecPlot(final_df))
  }

  rowData(exp) <- final_df
  if (filter) exp <- filterIsotopes(exp)
  return(exp)
}

#' @title Remove identified isotopes from the experiment
#' @param exp SummarizedExperiment with identified isotopes
#' @export
filterIsotopes <- function(exp){
  exp[-grep("isotope", rowData(exp)$isotopic_status), ]
}
