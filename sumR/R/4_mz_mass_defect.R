# *** MassDefect -----------------------------------------------------
#' MZ-OBJ MassDefect
#'
#' Mass Defect Calculation.\cr
#'     Zeros any intensities that are lower than MD-Calc.
#'
#' @param input_mz_obj \cr
#'   DataFrame : Input MZ-Obj
#' @param plot \cr
#'   Boolean   : To plot or not to plot.
#' @param magicNumber1 \cr
#'   Float     : MagicNumber1 in "MN1 * data + MN2"
#' @param magicNumber2 \cr
#'   Float     : MagicNumber2 in "MN1 * data + MN2"
#'
#' Dependencies : ggplot2, ggpubr, dplyr
#' @return Returns an MZ-OBJ
#' @export

mz_mass_defect <- function(input_mz_obj, plot = TRUE, magicNumber1 = 0.00112, magicNumber2 = 0.01953) {
  input_mz_obj$MD <-  as.numeric(input_mz_obj$mz) %% 1
  input_mz_obj$mz_filter <- magicNumber1 * as.numeric(input_mz_obj$mz) + magicNumber2

  #(revisit if original SUMR hmdb/McMillan is important)
  md_filtered <- input_mz_obj[which(input_mz_obj$MD <= input_mz_obj$mz_filter), ]

  mz_removed <- local.mz_log_removed_rows(input_mz_obj,md_filtered,
                  "sumR::mz_mass_defect")["mz_removed"]

  tryCatch({
    if(plot) {
      ggpubr::ggscatter(md_filtered, x="mz", y="MD",
                        fill     = "MD",
                        color    = "MD",
                        xlab     = "m/z",
                        ylab     = "MD",
                        title    = paste("MassDefect",local.mz_polarity_guesser(input_mz_obj),sep="-"),
                        subtitle = paste("Rows removed: ", mz_removed),
                        legend = "none"
      )

      local.save_plot(paste("MassDefect",local.mz_polarity_guesser(input_mz_obj),sep="-"))
    }
  }, error = function(e) {
      print("WARN: mz_mass_defect did not filter out anything to plot")
      print(e)
  }, finally = {
    return(dplyr::select(md_filtered, -MD, -mz_filter))
  })
}
