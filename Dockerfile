FROM  bioconductor/bioconductor_docker:RELEASE_3_17-R-4.3.0 AS bioc_base

RUN R -e 'BiocManager::install("mzR")'
RUN curl https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip -o hmdb_metabolites.zip
RUN unzip hmdb_metabolites.zip; mv hmdb_metabolites.xml /usr/hmdb_metabolites.xml

RUN curl "https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip" -o lipids.zip
RUN unzip lipids.zip; mv structures.sdf /usr/lipids.sdf


FROM bioc_base AS package_base

RUN R -e 'remotes::install_version("dplyr",        version = "1.1.2")'
RUN R -e 'remotes::install_version("ggplot2",      version = "3.4.2")'
RUN R -e 'remotes::install_version("tibble",       version = "3.2.1")'
RUN R -e 'remotes::install_version("tidyr",        version = "1.3.0")'
RUN R -e 'remotes::install_version("pbapply",      version = "1.7-2")'
RUN R -e 'remotes::install_version("devtools",     version = "2.4.5")'
RUN R -e 'remotes::install_version("pheatmap",     version = "1.0.12")'
RUN R -e 'remotes::install_version("XML",          version = "3.99-0.14")'
RUN R -e 'remotes::install_version("randomForest", version = "4.7-1.1")'
RUN R -e 'remotes::install_version("caret",        version = "6.0-94")'
RUN R -e 'remotes::install_version("caTools",      version = "1.18.2")'
RUN R -e 'remotes::install_version("AICcmodavg",   version = "2.3-2")'
RUN R -e 'remotes::install_version("ggpubr",       version = "0.6.0")'
RUN R -e 'remotes::install_version("pROC",         version = "1.18.5")'

#User Tools
RUN R -e 'remotes::install_version("readxl",   version = "1.4.3")'
RUN R -e 'remotes::install_version("testthat", version = "3.2.0")'


FROM package_base AS db_base
COPY scripts/db_massager.sh db_massager.sh
RUN sh db_massager.sh

FROM db_base
COPY MeDUSA MeDUSA
WORKDIR MeDUSA

## Create a more manageable hmdb xml file

RUN chown rstudio /home/rstudio/*.xml

RUN R -e 'devtools::document()'
RUN R -e 'devtools::install(dependencies="never")'


### TO BUILD
# docker build . -f Dockerfile-localR -t lacdr/medusa
### TO BUILD without MeDUSA (i.e. to develop it)
# docker build . -f Dockerfile-localR -t lacdr/medusa --target package_base

### TO RUN:
### TO RUN RSTUDIO: user=rstudio, pwd=medusa
# docker run -e PASSWORD=medusa -p 8787:8787 -v .:/home/rstudio/local lacdr/medusa
# Navigate to localhost:8787 in your browser

### TO RUN R terminal:
#[Unix/Mac] docker run --rm -it --name localR -v $(pwd):/local  lacdr/medusa /bin/bash
#[Win?]     docker run --rm -it --name localR -v .:/local  lacdr/medusa /bin/bash

### To enter a running container (use "docker ps" to make sure there is only one)
# docker exec --it localR /bin/bash
#[Unix/Mac] docker exec -it $(docker ps | grep medusa | awk '{print $1}') /bin/bash
