#!/bin/bash
# Product generation script for hdwx-modelplotter
# Created 15 September 2021 by Sam Gardner <stgardner4@tamu.edu>

rm plotcmds.txt
if [ -f lockfile ]
then
  echo "lockfile found, exiting"
  exit
fi
touch lockfile

#models=("hrrr", "nam", "gfs", "namnest")
~/mambaforge/envs/HDWX/bin/python3 modelFetch.py gfs
plotcmdStr=`cat plotcmds.txt`
IFS=$'\n' plotcmdArr=($plotcmdStr)
counter=0
for plotcmd in "${plotcmdArr[@]}"
do
    time eval "$plotcmd" &
    procpids[${counter}]=$!
    ((counter=counter+1))
    while [ ${#procpids[@]} == 12 ]
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
rm lockfile
