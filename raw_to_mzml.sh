#!bin/bash

# This script is used to convert .RAW DI-MS data into .mzML
# Please install Docker first and ensure the `dockerd` is running
# Call this script with the -d argument to specify a folder where .RAW files are located
# This script will place .mzML in that folder with the same name
# If no argument is given, the script will assume the .RAW is in the current directory from which this script is called

# Set base variables for pipeline
DATA_FOLDER="$PWD"

# Look into arguments given
while getopts d: flag 
do
    case "${flag}" in
        d) DATA_FOLDER=${OPTARG};;
    esac
done

docker run -it --rm -e WINEDEBUG=-all -v "$DATA_FOLDER":/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/*.RAW --simAsSpectra --32