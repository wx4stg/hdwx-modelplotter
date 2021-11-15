#!/bin/bash
# Product generation script for hdwx-modelplotter
# Created 15 September 2021 by Sam Gardner <stgardner4@tamu.edu>


rm -rf output/
rm -rf inputGemFiles/
rm -rf modelData/
if [ ! -f gempakToCF.jar ]
then
    echo "Fetching unidata gempakToCF.jar..."
    curl -L https://github.com/mjames-upc/edu.ucar.unidata.edex.plugin.gempak/blob/master/gempakToCF.jar?raw=true -o gempakToCF.jar
fi
mkdir inputGemFiles/
echo "Transferring data from LDM..."
filesToXfr=`find /coriolis-ldm/gempak/model/hrrr/ -mmin -5 -type f -exec ls {} +`
for file in $filesToXfr
do
    cp $file ./inputGemFiles
done
echo "Converting files..."
mkdir modelData/
for inputFile in inputGemFiles/*
do
    outputFileName=`basename $inputFile .gem`
    outputFileName="$outputFileName.nc"
    java -jar gempakToCF.jar $inputFile modelData/$outputFileName 
done
mkdir -p output/gisproducts/hrrr/sfcT/
mkdir -p output/products/hrrr/sfcTwindMSLP/


~/miniconda3/envs/HDWX/bin/python3 modelPlot.py
