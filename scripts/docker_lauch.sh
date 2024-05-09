#!/bin/bash

docker run -e PASSWORD=sumr -p 8787:8787 -v .:/home/rstudio/ lacdr/sumr-build
