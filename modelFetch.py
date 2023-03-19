#!/usr/bin/env python3
# Script for interfacing with NOAA Operational Model Archive and Distribution System (NOMADS)
# Created 17 October 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from datetime import datetime as dt, timedelta
from os import path, remove
import json
import requests
from pathlib import Path
import sys
from ecmwf.opendata import Client
from time import sleep

# modelFetch.py <"gfs"/"nam"/"namnest"/"hrrr">
modelName = sys.argv[1]
euroVarList = {
    "t2m.grib2" : ["2t"],
    "sfcwind.grib2" : ["10u", "10v"],
    "sp.grib2" : ["msl"],
    "sfccomposite" : list(),
    "winds.grib2" : ["u", "v"],
    "heights.grib2" : ["gh"],
    "temps.grib2" : ["t"],
    "rh.grib2" : ["r"],
    "500vort" : list(),
    "jetisotachs" : list(),
    "850temps" : list(),
    "700rh" : list(),
    "4pnl" : list()
}
ncepVarList = {
    "t2m.grib2" : "&lev_2_m_above_ground=on&var_TMP=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # 2m Temperature
    "sfcwind.grib2" : "&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # 10m u and v
    "sp.grib2" : "&lev_surface=on&var_HGT=on&var_PRES=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # surface pressure and orography
    "sfccomposite" : "",
    "winds.grib2" : "&lev_250_mb=on&lev_500_mb=on&lev_850_mb=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # winds at pressure surface
    "refc.grib2" : "&var_REFC=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # composite simulated reflectivity
    "udh.grib2" : "&lev_5000-2000_m_above_ground=on&var_MXUPHL=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # updraft helicity
    "refccomposite" : "",
    "refd.grib2" : "&lev_1000_m_above_ground=on&var_REFD=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # 1km AGL simulated reflectivity
    "refdcomposite" : "",
    "heights.grib2" : "&lev_1000_mb=on&lev_250_mb=on&lev_500_mb=on&lev_850_mb=on&var_HGT=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # geopotential height at pressure surface
    "temps.grib2" : "&lev_850_mb=on&var_TMP=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # temperature at pressure surface
    "rh.grib2" : "&lev_700_mb=on&var_RH=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", # rel humidity at pressure surface
    "500vort" : "",
    "jetisotachs" : "",
    "850temps" : "",
    "700rh" : "",
    "4pnl" : ""

}
basePath = path.dirname(path.abspath(__file__))
client = Client(source="ecmwf")
def writeToCmd(stringToWrite):
    while path.exists(path.join(basePath, "plotter-is-reading")):
        sleep(.01)
    with open(path.join(basePath, "plotcmds.txt"), "a") as cmdw:
        cmdw.write(stringToWrite)
        cmdw.close()

def deleteLockOnCrash(exceptionType, exceptionVal, traceback):
    remove(path.join(basePath, path.join(basePath, "downloaderlock-"+modelName+".txt")))

def writeToStatus(stringToWrite):
    print(stringToWrite)
    stringToWrite = stringToWrite+"\n"
    with open(path.join(basePath, "status.txt"), "a") as statw:
        statw.write(stringToWrite)
        statw.close()

def fetchEuroModel(initRun, fHour, outputDir):
    if modelName == "ecmwf-hres":
        if initRun.hour in [0, 12]:
            requestedType = "fc"
            requestedStream = "oper"
        else:
            requestedType = "fc"
            requestedStream = "scda"
    elif modelName == "ecmwf-ens":
        requestedStream = "enfo"
    for filename, reqVariables in euroVarList.items():
        resDT = ""
        if len(reqVariables) > 0:
            try:
                res = client.retrieve(
                    type=requestedType,
                    stream=requestedStream,
                    date=initRun,
                    step=fHour,
                    param=reqVariables,
                    target=path.join(outputDir, filename)
                )
                resDT = res.datetime
            except Exception as e:
                print(str(e))
                return False
            writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun.strftime("%Y%m%d%H%M")+" "+str(fHour)+" "+filename.replace(".grib2", "")+"\n")
        else:
            resDT = initRun
            writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun.strftime("%Y%m%d%H%M")+" "+str(fHour)+" "+filename.replace(".grib2", "")+"\n")
    if resDT == initRun:
        return True

def fetchNcepModel(initRun, fHour, outputDir, templateStr):
    requestedForecastHour = str(f'{fHour:02}')
    requestedForecastHourLong = str(f'{fHour:03}')
    for filename, reqVariable in ncepVarList.items():
        if filename == "udh.grib2":
            if modelName not in ["namnest", "hrrr"]:
                continue
        if "grib2" in filename:
            urlToFetch = templateStr.replace("<REQUESTED_VARIABLE>", reqVariable).replace("<MODEL_INIT_TIME>", initRun.strftime("%H")).replace("<MODEL_INIT_DATE>", initRun.strftime("%Y%m%d")).replace("<FHOUR_LONG>", requestedForecastHourLong).replace("<FHOUR_SHORT>", requestedForecastHour)
            if modelName == "hrrr" and reqVariable == "&lev_700_mb=on&var_RH=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F":
                urlToFetch = urlToFetch.replace("&lev_700_mb=on&var_RH=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F", "&lev_700_mb=on&var_DPT=on&var_TMP=on&subregion=&leftlon=-130&rightlon=-60&toplat=50&bottomlat=20&dir=%2F")
            modelData = requests.get(urlToFetch)
            if "GRIB" in modelData.text:
                with open(path.join(outputDir, filename), "wb") as f:
                    f.write(modelData.content)
                writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun.strftime("%Y%m%d%H%M")+" "+str(fHour)+" "+filename.replace(".grib2", "")+"\n")
            else:
                return False
        else:
            writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun.strftime("%Y%m%d%H%M")+" "+str(fHour)+" "+filename.replace(".grib2", "")+"\n")
    return True

if __name__ == "__main__":
    sys.excepthook = deleteLockOnCrash
    if path.exists(path.join(basePath, "downloaderlock-"+modelName+".txt")):
        writeToStatus("Downloader is locked for model "+modelName+", goodbye")
        exit()
    else:
        with open(path.join(basePath, "downloaderlock-"+modelName+".txt"), "a") as lockWrite:
            lockWrite.write(dt.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    writeToStatus("Statring download routine for model "+modelName)
    if path.exists(path.join(basePath, "firstPlotDT.txt")):
        readFirstPlotFile = open(path.join(basePath, "firstPlotDT.txt"), "r")
        firstPlotTime = dt.strptime(readFirstPlotFile.read(), "%Y%m%d%H%M")
        readFirstPlotFile.close()
    else:
        firstPlotTime = dt.utcnow()
        writeFirstPlotFile = open(path.join(basePath, "firstPlotDT.txt"), "w")
        writeFirstPlotFile.write(firstPlotTime.strftime("%Y%m%d%H%M"))
        writeFirstPlotFile.close()
    if modelName == "gfs":
        productTypeBase = 300
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 120, 1)) + list(range(120, 385, 3))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t<MODEL_INIT_TIME>z.pgrb2.0p25.f<FHOUR_LONG><REQUESTED_VARIABLE>gfs.<MODEL_INIT_DATE>%2F<MODEL_INIT_TIME>%2Fatmos"
    elif modelName == "nam":
        productTypeBase = 500
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 36, 1)) + list(range(36, 85, 3))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam.pl?file=nam.t<MODEL_INIT_TIME>z.awphys<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
    elif modelName == "namnest":
        productTypeBase = 600
        longRuns = list(range(0, 19, 6))
        shortRuns = list()
        fHoursLongRun = list(range(0, 61, 1))
        fHoursShortRun = list()
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl?file=nam.t<MODEL_INIT_TIME>z.conusnest.hiresf<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
    elif modelName == "hrrr":
        productTypeBase = 800
        longRuns = list(range(0, 19, 6))
        shortRuns = [hour for hour in list(range(0, 24, 1)) if hour not in longRuns]
        fHoursLongRun = list(range(0, 49, 1))
        fHoursShortRun = list(range(0, 19, 1))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?file=hrrr.t<MODEL_INIT_TIME>z.wrfsfcf<FHOUR_SHORT>.grib2<REQUESTED_VARIABLE>hrrr.<MODEL_INIT_DATE>%2Fconus"
    elif modelName == "ecmwf-hres":
        productTypeBase = 1000
        longRuns = list(range(0, 13, 12))
        shortRuns = [hour for hour in list(range(0, 19, 6)) if hour not in longRuns]
        fHoursLongRun = list(range(0, 144, 3)) + list(range(144, 241, 6))
        fHoursShortRun = list(range(0, 91, 3))
    today = dt.utcnow().replace(hour=0, minute=0, second=0, microsecond=0)
    yesterday = today - timedelta(days=1)
    runsAndFHours = dict()
    for dateInt in [int(yesterday.strftime("%Y%m%d0000")), int(today.strftime("%Y%m%d0000"))]:
        for longRun in longRuns:
            dtKey = dateInt + (longRun * 100)
            dtKeyDt = dt.strptime(str(dtKey), "%Y%m%d%H%M")
            if dtKeyDt < dt.utcnow() and dtKeyDt > firstPlotTime:
                runsAndFHours[dtKey] = fHoursLongRun
        for shortRun in shortRuns:
            dtKey = dateInt + (shortRun * 100)
            dtKeyDt = dt.strptime(str(dtKey), "%Y%m%d%H%M")
            if dtKeyDt < dt.utcnow() and dtKeyDt > firstPlotTime:
                runsAndFHours[dtKey] = fHoursShortRun
    for run in runsAndFHours.copy().keys():
        runfile = path.join(basePath, "output", "metadata", "products", str(productTypeBase), str(run)+".json")
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
    if path.exists(path.join(basePath, "plotcmds.txt")):
        with open(path.join(basePath, "plotcmds.txt"), "r") as pltcmdfile:
            queuedPlotCommands = pltcmdfile.read()
    else:
        queuedPlotCommands = ""
    for runInt, hours in runsAndFHours.items():
        failedFHoursThisRun = list()
        runDT = dt.strptime(str(runInt), "%Y%m%d%H%M")
        for fhour in hours:
            plotCmdString = sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+str(runInt)+" "+str(fhour)
            frameFetchDidSucceed = True
            if plotCmdString not in queuedPlotCommands:
                frameFetchDidSucceed = False
                writeToStatus("Fetching "+modelName+" init at "+runDT.strftime("%Y-%m-%d %H:%M")+" forecast hour "+str(fhour))
                outputDir = path.join(basePath, "modelData", modelName, runDT.strftime("%Y%m%d"), runDT.strftime("%H"), str(fhour))
                Path(outputDir).mkdir(parents=True, exist_ok=True)
                if "ecmwf" in modelName:
                    frameFetchDidSucceed = fetchEuroModel(runDT, fhour, outputDir)
                else:
                    frameFetchDidSucceed = fetchNcepModel(runDT, fhour, outputDir, templateString)
            if frameFetchDidSucceed == False:
                failedFHoursThisRun.append(fhour)
            if len(failedFHoursThisRun) > 3:
                break
    remove(path.join(basePath, "downloaderlock-"+modelName+".txt"))
    print("Finished downloading for model "+modelName)