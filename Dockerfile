### This image is quite heavy. Can take over an hour to build.
FROM  bioconductor/bioconductor_docker:RELEASE_3_17-R-4.3.0 AS bioc_base

RUN R -e 'BiocManager::install("mzR")'
COPY scripts/db_download.sh scripts/db_massager.sh ./
RUN sh /db_download.sh; sh /db_massager.sh

RUN R -e 'remotes::install_version("dplyr",        version = "1.1.2")'   \
&&  R -e 'remotes::install_version("ggplot2",      version = "3.4.2")'   \
&&  R -e 'remotes::install_version("tibble",       version = "3.2.1")'   \
&&  R -e 'remotes::install_version("tidyr",        version = "1.3.0")'   \
&&  R -e 'remotes::install_version("pbapply",      version = "1.7-2")'   \
&&  R -e 'remotes::install_version("devtools",     version = "2.4.5")'   \
&&  R -e 'remotes::install_version("pheatmap",     version = "1.0.12")'  \
&&  R -e 'remotes::install_version("XML",          version = "3.99-0.14")' \
&&  R -e 'remotes::install_version("randomForest", version = "4.7-1.1")' \
&&  R -e 'remotes::install_version("caret",        version = "6.0-94")'  \
&&  R -e 'remotes::install_version("caTools",      version = "1.18.2")'  \
&&  R -e 'remotes::install_version("AICcmodavg",   version = "2.3-2")'   \
&&  R -e 'remotes::install_version("Matrix",       version = "1.6.0")'   \
&&  R -e 'remotes::install_version("ggpubr",       version = "0.6.0")'   \
&&  R -e 'remotes::install_version("pROC",         version = "1.18.5")'  \
&&  R -e 'remotes::install_version("readxl",       version = "1.4.3")'   \
&&  R -e 'remotes::install_version("testthat",     version = "3.2.1.1")'

FROM bioc_base AS medusa_test
COPY MeDUSA MeDUSA
WORKDIR MeDUSA
RUN chown rstudio /home/rstudio/*.xml
RUN R -e 'devtools::document()' \
&&  R -e 'devtools::test()' \
&&  R -e 'devtools::install(dependencies="never")'

FROM bioc_base AS medusa_image
COPY MeDUSA MeDUSA
WORKDIR MeDUSA
RUN chown rstudio /home/rstudio/*.xml

RUN R -e 'devtools::document()' \
&&  R -e 'devtools::install(dependencies="never")'

### TO BUILD
# docker build . -f Dockerfile -t lacdr/medusa

### TO RUN RSTUDIO: user=rstudio, pwd=medusa
# docker run -e PASSWORD=medusa -p 8787:8787 -v .:/home/rstudio/local lacdr/medusa
# Navigate to localhost:8787 in your browser
# Note, this creates a shared volume. So best to run this from the directory of your data.
#     (But not so high of a directory that it will consume unnecessary resources)

### TO RUN R terminal:
#[Unix/Mac] docker run --rm -it --name localR -v $(pwd):/local  lacdr/medusa /bin/bash
#[Win?]     docker run --rm -it --name localR -v .:/local  lacdr/medusa /bin/bash
### To enter a running container (use "docker ps" to make sure there is only one)
# docker exec --it localR /bin/bash
#[Unix/Mac] docker exec -it $(docker ps | grep medusa | awk '{print $1}') /bin/bash
