#!/bin/bash
# Product generation script for hdwx-modelplotter
# Created 15 September 2021 by Sam Gardner <stgardner4@tamu.edu>


if [ ! -f gempakToCF.jar ]
then
    echo "Fetching unidata gempakToCF.jar..."
    curl -L https://github.com/mjames-upc/edu.ucar.unidata.edex.plugin.gempak/blob/master/gempakToCF.jar?raw=true -o gempakToCF.jar
fi
echo "Converting files..."
for inputFile in inputGemFiles/*
do
    outputFileName=`basename $inputFile .gem`
    outputFileName="$outputFileName.nc"
    java -jar gempakToCF.jar $inputFile $outputFileName 
done
mkdir -p output/gis/
