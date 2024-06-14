#!/bin/bash

docker run -e PASSWORD=medusa -p 8787:8787 -v .:/home/rstudio/ lacdr/medusa
