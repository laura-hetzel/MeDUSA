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


# Still in development
#' @title Annotation for adducts incl isotopes
#' @description A list of annotation information of isotopes
#' @usage annotate_iso(mz_vector, ppm , z )
#' @param mz_vector mz value vector
#' @param ppm error part per million
#' @param z charge values
annotate_add <- function(mz_vector){
  x <- massdiff(mz_vector)
  y <- adductMatch(x, add = mass2adduct::adducts2, ppm = 5)
}

#' @title Annotation for isotopes
#' @description A list of annotation information of isotopes
#' @usage annotate_iso(mz_vector, ppm , z )
#' @param mz_vector mz value vector
#' @param ppm error part per million
#' @param z charge values
#' @importFrom Rdisop decomposeMass
annotate_iso <- function(mz_vector, ppm , z ){
  list_annotate <- c()
  for (i in 1:length(mz_vector)) {
    list_annotate[[i]] <- Rdisop::decomposeMass(mz_vector[i], ppm = ppm,z = z, maxisotopes = 10, minElements = "C0", maxElements = "C999999")
  }
  return(list_annotate)
}


#' @title Mass spectrometer visualization
#' @description Gives a graphical representation of the mono-isotopic and isotopic peaks
#' @usage mass_spec_plot(final_df_vis)
#' @param final_df_vis dataframe of the mono-isotopic and isotopic peaks pair's with id column
#' @importFrom  plotly plot_ly highlight_key add_markers add_segments  layout rangeslider ggplotly highlight
#' @importFrom dplyr %>%
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


# Still in development
# Maybe replaced later
#' @title Imputation method for missing intensity values
#' @description Imputation method for missing intensity values
impute.KNN.obs.sel <- function(dat, # incomplete data matrix
                               cor.var.sel = 0.2, # correlation threshold for variable pre-selection
                               K=5, # number of neighbors
                               verbose=T) {

  datimp <- dat

  # incomplete observations
  incom.obs <- which(apply(dat,1,function(x) any(is.na(x))))
  if (verbose) message(paste0("Number of imcomplete observations: ", length(incom.obs)))

  # incomplete variables
  incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
  if (verbose) message("Proceeding obs... ")

  # Pearson correlation for variable pre-selection
  Cor <- cor(dat,use = "p")

  # Calculate distance
  D2list <- lapply(incom.vars, function(j) {
    varsel <- which(abs(Cor[j,]) > cor.var.sel)
    if (length(varsel) > 10) varsel <- order(abs(Cor[j,]),decreasing = T)[1:11]
    if (length(varsel) < 5) varsel <- order(abs(Cor[j,]),decreasing = T)[1:6]
    D2 <- as.matrix(dist(scale(dat[,varsel])),upper = T,diag = T)
    if (any(is.na(D2))) {
      D2a <- as.matrix(dist(scale(dat)),upper = T,diag = T)*sqrt(length(varsel)/ncol(dat))
      D2[is.na(D2)] <- D2a[is.na(D2)] }
    diag(D2) <- NA
    D2})
  names(D2list) <- incom.vars

  # get neighbors and impute by their weighted average
  for (i in incom.obs) {
    if (verbose) cat(paste(i," "))
    comvars <-  complete.cases(as.numeric(dat[i,]))
    dattmp <-  dat
    for (j in which(!comvars)) {
      D2 <- D2list[[as.character(j)]]
      if (any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last = NA)
        KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(dat[ii,j])))]
      }  else KNNids  <- NULL

      if (!is.null(KNNids)) KNNids_sel <- intersect(KNNids[1:min(K,length(KNNids))],KNNids_naomit)
      if (length(KNNids_sel) < 1) KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))] else
        if (length(which(sapply(KNNids_sel,function(ii) !is.na(dat[ii,j])))) < floor(K/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(K/2),length(KNNids_naomit))]
        if (any(!is.na(D2[i,])) & length(KNNids)>= 1) {
          dattmp_sel <- dattmp[KNNids_sel,j]
          dattmp[i,j] <- sum(dattmp_sel*exp(-D2[i,KNNids_sel]),na.rm = T)/sum(exp(-D2[i,KNNids_sel])[!is.na(dattmp_sel)],na.rm=T) }
    }
    datimp[i,] <- dattmp[i,]
  }

  # if there are still NAs, take mean
  datimp <- apply(datimp,2, function(x) {
    if (any(is.na(x))) x[is.na(x)] <- mean(x,na.rm = T)
    x})



  datimp}


#---------------------------------------------------------------------------------------------------
#' @title Isotope tagging
#' @description Tags isotopes based on intensity ratio and difference in mass
#' @usage Isotope_tagging(df, iso_diff , ppm , by_step )
#' @param df Dataframe must contain columns c("mz", "intensity") with "mz" as first column
#' @param ppm Parts Per Million Tolerance
#' @param z charge
#' @importFrom  data.table fintersect
#' @importFrom  dplyr mutate tibble filter select mutate rename
#' @importFrom tibble tibble
#' @importFrom tidyselect everything
#' @export
#'
isotope_tagging <- function(df , iso_diff_da, ppm , by_step, z){
  #-----------------------------------------------------------------------------------
  ## Start body part of function BAIT

  #keep orginal copy
  df1 <- df

  #creating new data frame with mean intensity column # note now complete case but will be replaced by imputation method
  df <- data.frame(mz = df[,1], intensity = rowMeans(df[,-1], na.rm = TRUE))
  #setting mz & intensity to vectors
  mz_vector <-  df$mz

  #creating new df with mass_error , lower bound and upperbound
  df <- mutate(df,
               mass_error = ppm_to_dalton(mz, ppm),
               min_error  = mz - mass_error,
               max_error  = mz + mass_error)


  #Annotation
  adducts <- annotate_add(mz_vector)


#Isotope identification
  # (1) modifying the df
  add_a <- filter(df, df$mz %in% adducts$A)
  add_a <- add_a[,1:2]
  add_a$intensity_A <- add_a$intensity
  add_a$intensity <- NULL

  adduct_A <- merge(x = adducts, y = add_a, by.x = "A", by.y = "mz", all.x = TRUE)

  add_b <- filter(df, df$mz %in% adducts$B)
  add_b <- add_b[,1:2]
  add_b$intensity_B <- add_b$intensity
  add_b$intensity <- NULL

  adduct <- merge(x = adducts, y = add_b, by.x = "B", by.y = "mz", all.x = TRUE)

  adduct$intensity_A <- adduct_A$intensity_A
  # (2) filter for isotope elements
  adducts_isotope <- filter(adduct, adduct$matches %in% c("13C"))

  # (3) moleculist_isotopes_list
  isotope_molecules <- NULL
  length_iso <- nrow(adducts_isotope)
  for (i in 1:length_iso) {
    masses <- c(adduct$A[i], adduct$B[i])
    intensities <- c(100, ((adduct$intensity_B[i]/adduct$intensity_A[i])*100) )
    isotope_molecules[[i]] <- decomposeIsotopes(masses , intensities,  ppm = 5 , z = 1)
  }

  # (4) filtering for validity
  isotope_valid_list <- NULL
  isotope_molecules[sapply(isotope_molecules, is.null)] <- NULL
  length_iso_v <- length(isotope_molecules)
  for (i in 1:length_iso_v) {
    if (isotope_molecules[[i]]$valid == "Valid") {
      isotope_valid_list[[i]] <- lapply(isotope_molecules[[i]], `[[`, 1)
    }
  }

  # (5) Creating data fram from list
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

  # (6) Filter for formula with "C" in it
  isotope_valid_df <- isotope_valid_df[grep("^C",isotope_valid_df$formula),]

  # (7) Change to numeric
  isotope_valid_df$exactmass <- as.numeric(isotope_valid_df$exactmass)
  isotope_valid_df$isotope_mass <- as.numeric(isotope_valid_df$isotope_mass)

  # (8) Merge together
  merge_df <- sqldf("select * from isotope_valid_df left join df
      on (isotope_valid_df.isotope_mass > df.min_error and isotope_valid_df.isotope_mass <= df.max_error) ")
  merge_df <- na.omit(merge_df)
  merge_dff <- sqldf("select * from isotope_valid_df inner join df on (isotope_valid_df.exactmass > df.min_error and isotope_valid_df.exactmass <= df.max_error) ")
  merge_dff <- na.omit(merge_dff)

  final <- merge(x = merge_df, y = merge_dff, by = "formula", all = TRUE)
  final <- na.omit(final)
  #(9)Isotopic Status
  final_df <- select(mutate(df,
                            Isotopic_status = case_when(
                              mz_vector  %in% final$mz.x~ "c13_isotope",
                              mz_vector  %in% final$mz.y ~ "molecular_ion")),
                     Isotopic_status ,
                     everything())

  #Adducts
  adduct_b <- adduct[,c(1,5)]
  final_df <- merge(final_df,adduct_b, by.x = "mz", by.y = "B", all = TRUE)

  #plotting
  print(mass_spec_plot(final_df_vis))


  #adding colums to original dataframe
  df1$Isotopic_Status <- final_df$Isotopic_status
  df1$min_error <- df$min_error
  df1$max_error <- df$max_error
  df1$Adduct <- final_df$matches




  return(df1)
}






