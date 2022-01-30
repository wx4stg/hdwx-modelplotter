#!/bin/bash
# HDWX operational submodule resetter script
# Created 10 Janurary 2022 by Sam Gardner <stgardner4@tamu.edu>

rm -rf output/
rm -rf radarData/
rm downloaderlock-*
rm plotter-is-reading
rm plotterlock
rm status.txt
rm plotcmds.txt
rm firstPlotDT.txt