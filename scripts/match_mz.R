#' Title
#'
#' @param df_raw
#' @param ppm
#' @param ionmode
#' @param base.dbname name of the database that should be used, e.g. hmdb, lipidmaps,
#'
#' @param databaseDir path to the folder containing the downloaded database.
#' See example on how to download one.
#'
#' @return
#' @export
#'
#' @importFrom MetaDBparse searchMZ
#' @importFrom tidyr separate
#'
#' @examples
#' \dontrun{
#' ## instead of a temporary dir, define a fixed path
#' database_dir = tempdir()
#' dbname = "lipidmaps"
#' buildBaseDB(outfolder = database_dir, dbname = dbname, test = FALSE)
#'
#' matches <- match_mz_with_database(df_raw = df_raw, database_dir = database_dir,
#' base.dbname = dbname, ppm = 5, ionmode = c("positive"))
#' }
match_mz_with_database <- function(df_raw, database_dir, base.dbname = "hmdb",
                                   ppm = 5, ionmode = c("positive")){

    ## forward searching
    df_sep <- separate(df_raw, col = mz, into = c("ion_mode", "mz"), sep = "_")
    df_sep2 <- df_sep %>%
        mutate(mz == as.numeric(mz)) %>%
        # I guess pm is short for positive mode?
        # If so need to a condition here that switches between positive and
        # negative mode
        filter(ion_mode == "pm")

    matches <- searchMZ(mzs = df_sep2$mz,
                        ionmodes = ionmode,
                        outfolder = database_dir,
                        base.dbname = base.dbname,
                        ppm = ppm,
                        ext.dbname = "extended")

    return(matches)
}
