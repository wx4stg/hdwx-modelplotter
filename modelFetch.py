#!/usr/bin/env python3
# Script for interfacing with NOAA Operational Model Archive and Distribution System (NOMADS)
# Created 17 October 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from datetime import datetime as dt, timedelta, UTC
from os import path, listdir
import json
import requests
import sys
from ecmwf.opendata import Client
import tempfile
from bs4 import BeautifulSoup
import numpy as np
from modelPlot import plot_all


euroVarList = ["2t", "2d", "10u", "10v", "msl", "u", "v", "gh", "t", "r"]

ncepVarListMaster = (
    # levels
    "&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on&lev_5000-2000_m_above_ground=on&lev_1000_m_above_ground=on" # AGL
    "&lev_250_mb=on&lev_500_mb=on&lev_700_mb=on&lev_850_mb=on&lev_1000_mb=on&lev_entire_atmosphere_(considered_as_a_single_layer)=on" # Pressure
    # variables
    "&var_TMP=on" # Temperature
    "&var_DPT=on" # dewpoint
    "&var_RH=on" # rel humidity
    "&var_UGRD=on&var_VGRD=on" # 10m u and v
    "&var_PRES=on" # pressure
    "&var_HGT=on" # geopotential height + sfc orography
    "&var_REFC=on" # composite simulated reflectivity
    "&var_REFD=on" # simulated reflectivity
    "&var_MXUPHL=on" # updraft helicity
    # subregion
    "&subregion=&leftlon=-144.5&rightlon=-44.5&toplat=54.5&bottomlat=14.5&dir=%2F"
)

modelAvailableChecks = {
    'gfs' : {
        'url' : 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.%Y%m%d/%H/atmos/',
        'prefix' : 'gfs.t%Hz.pgrb2.0p25.f',
        'count' : 3
    },
    'nam' : {
        'url' : 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.%Y%m%d/',
        'prefix' : 'nam.t%Hz.awphys',
        'count' : 2
    },
    'namnest' : {
        'url' : 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.%Y%m%d/',
        'prefix' : 'nam.t%Hz.conusnest.hiresf',
        'count' : 2
    },
    'hrrr' : {
        'url' : 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/hrrr.%Y%m%d/conus/',
        'prefix' : 'hrrr.t%Hz.wrfsfcf',
        'count' : 2
    },
    'ecmwf-hres' : {
        'url' : 'https://data.ecmwf.int/forecasts/%Y%m%d/%Hz/ifs/0p25/oper/',
        'prefix' : '%Y%m%d%H0000-',
        'split' : 1
    }
}


basePath = path.dirname(path.abspath(__file__))

def fetchEuroModel(initRun, fHour):
    global euroVarList
    ecmwf = Client(source="ecmwf")
    if initRun.hour in [0, 12]:
        requestedType = "fc"
        requestedStream = "oper"
    else:
        requestedType = "fc"
        requestedStream = "scda"
    should_gis = ("--no-gis" not in sys.argv)
    with tempfile.NamedTemporaryFile(delete=False) as f:
        ecmwf.retrieve(type=requestedType, stream=requestedStream, date=initRun,
            step=fHour, resol="0p25",param=euroVarList,target=f.name)
        plot_all([f.name], 'ecmwf-hres', should_gis)

def fetchNcepModel(initRun, fHour, templateStr):
    global ncepVarListMaster
    ncepVarList = ncepVarListMaster
    requestedForecastHour = str(f'{fHour:02}')
    requestedForecastHourLong = str(f'{fHour:03}')
    modelName = None
    if 'gfs' in templateStr:
        modelName = 'gfs'
    elif 'nest' in templateStr:
        modelName = 'namnest'
    elif 'nam' in templateStr:
        modelName = 'nam'
    elif 'hrrr' in templateStr:
        modelName = 'hrrr'
    if 'nest' not in templateStr and 'hrrr' not in templateStr:
        ncepVarList = ncepVarList.replace("&lev_5000-2000_m_above_ground=on&lev_1000_m_above_ground=on", "").replace("&var_MXUPHL=on", "")
    if 'nam' not in templateStr:
        ncepVarList = ncepVarList.replace("_(considered_as_a_single_layer)", "")
    if 'hrrr' in templateStr or 'nest' in templateStr:
        ncepVarList = ncepVarList.replace("&lev_250_mb=on", "")
        hrrr_upper_air_template = templateStr.replace("<REQUESTED_VARIABLE>", "&var_UGRD=on&var_VGRD=on&lev_250_mb=on&subregion=&leftlon=-144.5&rightlon=-44.5&toplat=54.5&bottomlat=14.5&dir=%2F").replace("<MODEL_INIT_TIME>", initRun.strftime("%H")).replace("<MODEL_INIT_DATE>", initRun.strftime("%Y%m%d")).replace("<FHOUR_LONG>", requestedForecastHourLong).replace("<FHOUR_SHORT>", requestedForecastHour)
        hrrr_upper_air_data = requests.get(hrrr_upper_air_template)
        if hrrr_upper_air_data.text.startswith('GRIB'):
            with tempfile.NamedTemporaryFile(delete=False) as f:
                f.write(hrrr_upper_air_data.content)
                hrrr_upper_air_data = [f.name]
        else:
            if 'Data file is not present:' in hrrr_upper_air_data.text:
                return
            raise ValueError(f'Invalid HRRR upper-air data: {hrrr_upper_air_data.text}')
    else:
        hrrr_upper_air_data = []
    urlToFetch = templateStr.replace("<REQUESTED_VARIABLE>", ncepVarList).replace("<MODEL_INIT_TIME>", initRun.strftime("%H")).replace("<MODEL_INIT_DATE>", initRun.strftime("%Y%m%d")).replace("<FHOUR_LONG>", requestedForecastHourLong).replace("<FHOUR_SHORT>", requestedForecastHour)
    modelData = requests.get(urlToFetch)
    if modelData.text.startswith('GRIB'):
        should_gis = ("--no-gis" not in sys.argv)
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(modelData.content)
            datasets = [f.name]
            datasets.extend(hrrr_upper_air_data)
            plot_all(datasets, modelName, should_gis)
    else:
        if 'Data file is not present:' in modelData.text:
            return
        raise ValueError(f'Invalid model data: {modelData.text}')

def fetchFramesForModel(modelName):
    print("Statring download routine for model "+modelName)
    today = dt.now(UTC).replace(tzinfo=None)
    yesterday = today - timedelta(days=1)
    todayInt = int(today.replace(hour=0, minute=0, second=0, microsecond=0).strftime("%Y%m%d0000"))
    yesterdayInt = int(yesterday.replace(hour=0, minute=0, second=0, microsecond=0).strftime("%Y%m%d0000"))
    if modelName == "gfs":
        productTypeBase = 300
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 385, 3))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t<MODEL_INIT_TIME>z.pgrb2.0p25.f<FHOUR_LONG><REQUESTED_VARIABLE>gfs.<MODEL_INIT_DATE>%2F<MODEL_INIT_TIME>%2Fatmos"
        threshold = today - timedelta(hours=7)
    elif modelName == "nam":
        productTypeBase = 500
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 36, 1)) + list(range(36, 85, 3))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam.pl?file=nam.t<MODEL_INIT_TIME>z.awphys<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
        threshold = today - timedelta(hours=7)
    elif modelName == "namnest":
        productTypeBase = 600
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 61, 1))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl?file=nam.t<MODEL_INIT_TIME>z.conusnest.hiresf<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
        threshold = today - timedelta(hours=7)
    elif modelName == "hrrr":
        productTypeBase = 800
        longRuns = list(range(0, 19, 6))
        shortRuns = [hour for hour in list(range(0, 24, 1)) if hour not in longRuns]
        fHoursLongRun = list(range(0, 49, 1))
        fHoursShortRun = list(range(0, 19, 1))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?file=hrrr.t<MODEL_INIT_TIME>z.wrfsfcf<FHOUR_SHORT>.grib2<REQUESTED_VARIABLE>hrrr.<MODEL_INIT_DATE>%2Fconus"
        threshold = today - timedelta(hours=2)
    elif modelName == "ecmwf-hres":
        productTypeBase = 1000
        longRuns = list(range(0, 13, 12))
        shortRuns = [hour for hour in list(range(0, 19, 6)) if hour not in longRuns]
        fHoursLongRun = list(range(0, 144, 3)) + list(range(144, 241, 6))
        fHoursShortRun = list(range(0, 91, 3))
        threshold = today - timedelta(hours=14)
    runsAndFHours = dict()
    for dateInt in [yesterdayInt, todayInt]:
        for longRun in longRuns:
            dtKey = dateInt + (longRun * 100)
            dtKeyDt = dt.strptime(str(dtKey), "%Y%m%d%H%M")
            if dtKeyDt < today:
                if dtKeyDt > threshold:
                    runsAndFHours[dtKey] = fHoursLongRun
        for shortRun in shortRuns:
            dtKey = dateInt + (shortRun * 100)
            dtKeyDt = dt.strptime(str(dtKey), "%Y%m%d%H%M")
            if dtKeyDt < today:
                if dtKeyDt > threshold:
                    runsAndFHours[dtKey] = fHoursShortRun
    for run in runsAndFHours.copy().keys():
        if path.exists(path.join(basePath, "output", "metadata", "products")):
            runFileDirs = [d for d in listdir(path.join(basePath, "output", "metadata", "products")) if path.isdir(path.join(basePath, "output", "metadata", "products", d)) and int(d) <= productTypeBase+99 and int(d) >= productTypeBase]
            runFileDirs = sorted(runFileDirs, reverse=True)
            if len(runFileDirs) > 0:
                runfile = path.join(basePath, "output", "metadata", "products", runFileDirs[0], str(run)+".json")
                if path.exists(runfile):
                    with open(runfile) as jsonRead:
                        runMetadata = json.load(jsonRead)
                    if runMetadata["availableFrameCount"] == runMetadata["totalFrameCount"]:
                        runsAndFHours.pop(run)
                    else:
                        frmsToDelete = [frame["fhour"] for frame in runMetadata["productFrames"]]
                        oldArr = runsAndFHours[run]
                        newArr = [oldRun for oldRun in oldArr if oldRun not in frmsToDelete]
                        runsAndFHours[run] = newArr
        if len(runsAndFHours[run]) > 0:
            runDt = dt.strptime(str(run), "%Y%m%d%H%M")
            framesAvailable = []
            url = modelAvailableChecks[modelName]['url']
            if 'ecmwf' in modelName and runDt.hour % 12 != 0:
                url = url.replace('oper', 'scda')
            res = requests.get(runDt.strftime(url))
            soup = BeautifulSoup(res.text, features='lxml')
            this_prefix = runDt.strftime(modelAvailableChecks[modelName]['prefix'])
            for link in soup.find_all('a'):
                link = link.string
                if link.startswith(this_prefix):
                    if 'split' in modelAvailableChecks[modelName]:
                        framesAvailable.append(int(link.split('-')[modelAvailableChecks[modelName]['split']].replace('h', '')))
                    elif 'count' in modelAvailableChecks[modelName]:
                        framesAvailable.append(int(link.replace(this_prefix, '')[0:modelAvailableChecks[modelName]['count']]))
            runsAndFHours[run] = [frame for frame in runsAndFHours[run] if frame in np.unique(framesAvailable)]
    for runInt, hours in runsAndFHours.items():
        runDt = dt.strptime(str(runInt), "%Y%m%d%H%M")
        for fhour in hours:
            print("Fetching "+modelName+" init at "+runDt.strftime("%Y-%m-%d %H:%M")+" forecast hour "+str(fhour))
            if "ecmwf" in modelName:
                fetchEuroModel(runDt, fhour)
            else:
                fetchNcepModel(runDt, fhour, templateString)
    print("Finished processing for model "+modelName)

if __name__ == '__main__':
    for model in modelAvailableChecks.keys():
        fetchFramesForModel(model)
