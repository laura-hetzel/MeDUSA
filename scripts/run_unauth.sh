#!/usr/bin/with-contenv bash
## load /etc/environment vars first:

## Used to avoid auth-ing in browser.
## Pinned FROM  bioconductor/bioconductor_docker:RELEASE_3_17-R-4.3.0
## Future updates from the bioc built rocker image may break this flow

for line in $( cat /etc/environment ) ; do export $line > /dev/null; done
export USER=rstudio
exec /usr/lib/rstudio-server/bin/rserver --server-daemonize=0 --auth-none=1
