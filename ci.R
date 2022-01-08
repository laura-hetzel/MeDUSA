setwd("sumR")

install_if_needed <- function(package_to_install){
  if (!package_to_install %in% installed.packages()) {
    install.packages(package_to_install)
    print(sprintf("Installed %s", package_to_install))
  } else {
    print(sprintf("Found Installation of %s", package_to_install))
  }
}

ci_setup <- function(){
  options(repos = structure(BiocManager::repositories()))
  install_if_needed("devtools")
  install_if_needed("packrat")
  packrat::init(restart = FALSE,
                options = list(
                  external.packages = c(
                    "devtools", "roxygen2",
                    "remotes", "digest"
                  )
                ))
  packrat::restore()
}

ci_check <- function(){
  if (length(list.files(path = "R") > 0)) {
    ci_setup()
    devtools::check(error_on = "error")
  }
}

ci_coverage <- function(){
  if (length(list.files(path = "R") > 0)) {
    ci_setup()
    install_if_needed("covr")
    covr::package_coverage(type = c("tests", "examples"))
  }
}

