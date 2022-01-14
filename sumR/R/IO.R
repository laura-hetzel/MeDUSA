#'@title Extract identified features with intensities
#'@param data readMSData object
#'@param output Directory of where the data is stored
#'@importFrom xcms featureValues fillChromPeaks FillChromPeaksParam
#'featureDefinitions
#'@importFrom utils write.csv
extract_features <- function(data, output){
  write.csv(featureValues(data, value = "into"),
            file = sprintf("%s/Filled_Feature_values.csv", output))

  write.csv(featureDefinitions(data),
            file = sprintf("%s/Filled_Feature_definitions.csv", output))

  return(data)
}
