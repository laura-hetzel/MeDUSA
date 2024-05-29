#!/usr/env bash

R -e 'remotes::install_version("dplyr",        version = "1.1.2")'
R -e 'remotes::install_version("ggplot2",      version = "3.4.2")'
R -e 'remotes::install_version("ggplot2",      version = "3.4.2")'
R -e 'remotes::install_version("tibble",       version = "3.2.1")'
R -e 'remotes::install_version("tidyr",        version = "1.3.0")'
R -e 'remotes::install_version("pbapply",      version = "1.7-2")'
R -e 'remotes::install_version("devtools",     version = "2.4.5")'
R -e 'remotes::install_version("pheatmap",     version = "1.0.12")'
R -e 'remotes::install_version("XML",          version = "3.99-0.14")'
R -e 'remotes::install_version("randomForest", version = "4.7-1.1")'
R -e 'remotes::install_version("caret",        version = "6.0-94")'
R -e 'remotes::install_version("caTools",      version = "1.18.2")'
R -e 'remotes::install_version("AICcmodavg",   version = "2.3-2")'
R -e 'remotes::install_version("ggpubr",       version = "0.6.0")'
R -e 'remotes::install_version("pROC",         version = "1.18.5")'

#UserTools
R -e 'remotes::install_version("readxl",   version = "1.4.3")'
R -e 'remotes::install_version("testthat", version = "3.2.0")' ;
