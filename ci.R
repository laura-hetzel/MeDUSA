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
  if (!dir.exists("cache")) {
    dir.create("cache", showWarnings = F)
    print(sprintf("Creating cache folder at %s/cache", getwd()))
  }
  folder <- sprintf("%s/cache", getwd())
  .libPaths(folder)
  options(repos = structure(BiocManager::repositories()))
  install_if_needed("devtools")
  usethis::use_build_ignore("cache")
  devtools::install(upgrade = F)
}

ci_check <- function(){
  folder <- sprintf("%s/cache", getwd())
  .libPaths(folder)
  if (length(list.files(path = "R") > 0)) {
    devtools::check(error_on = "error")
  }
}

ci_coverage <- function(){
  folder <- sprintf("%s/cache", getwd())
  .libPaths(folder)
  if (length(list.files(path = "R") > 0)) {
    install_if_needed("covr")
    covr::package_coverage(type = c("tests", "examples"))
  }
}

