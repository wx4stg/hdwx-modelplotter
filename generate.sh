#!/bin/bash
# Product generation script for hdwx-modelplotter
# Created 15 September 2021 by Sam Gardner <stgardner4@tamu.edu>

if [ -f ../config.txt ]
then
    source ../config.txt
else
    condaEnvName="HDWX"
fi

if [ ! -d output/ ]
then
    mkdir output/
fi

if [ -f ~/mambaforge/envs/$condaEnvName/bin/python3 ]
then
    ~/mambaforge/envs/$condaEnvName/bin/python3 modelFetch.py ecmwf-hres &
    ~/mambaforge/envs/$condaEnvName/bin/python3 modelFetch.py namnest &
    ~/mambaforge/envs/$condaEnvName/bin/python3 modelFetch.py nam &
    ~/mambaforge/envs/$condaEnvName/bin/python3 modelFetch.py gfs &
    ~/mambaforge/envs/$condaEnvName/bin/python3 modelFetch.py hrrr &
fi
if [ -f ~/miniconda3/envs/$condaEnvName/bin/python3 ]
then
    ~/miniconda3/envs/$condaEnvName/bin/python3 modelFetch.py ecmwf-hres &
    ~/miniconda3/envs/$condaEnvName/bin/python3 modelFetch.py namnest &
    ~/miniconda3/envs/$condaEnvName/bin/python3 modelFetch.py nam &
    ~/miniconda3/envs/$condaEnvName/bin/python3 modelFetch.py gfs &
    ~/miniconda3/envs/$condaEnvName/bin/python3 modelFetch.py hrrr &
fi
counter=0
if [ ! -f plotterlock ]
then
    touch plotterlock
    while [ -f plotcmds.txt ]
    do
        touch plotter-is-reading
        read -r plotcmd < plotcmds.txt
        tail -n +2 plotcmds.txt > plotcmds.tmp && rm plotcmds.txt && mv plotcmds.tmp plotcmds.txt
        if ! grep -q '[^[:space:]]' "plotcmds.txt"
        then
            rm plotcmds.txt
        fi
        rm plotter-is-reading
        eval "$plotcmd" &
        procpids[${counter}]=$!
        ((counter=counter+1))
        while [ ${#procpids[@]} == 6 ]
        do
            for procpid in ${procpids[*]}
            do
                if ! kill -0 $procpid 2>/dev/null
                then 
                    procpids=(${procpids[@]/$procpid})
                fi
            done
        done
    done
    for procpid in ${procpids[*]}
    do
        wait $procpid
    done
    rm plotterlock
fi
if [ -f ~/mambaforge/envs/$condaEnvName/bin/python3 ]
then
    ~/mambaforge/envs/$condaEnvName/bin/python3 cleanup.py
fi
if [ -f ~/miniconda3/envs/$condaEnvName/bin/python3 ]
then
    ~/miniconda3/envs/$condaEnvName/bin/python3 cleanup.py
fi
