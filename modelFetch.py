#!/usr/bin/env python3
# Script for interfacing with NOAA Operational Model Archive and Distribution System (NOMADS)
# Created 17 October 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from datetime import datetime as dt, timedelta
from os import path, listdir
import json
import requests
from pathlib import Path
import subprocess
import sys

# modelFetch.py <"gfs"/"nam"/"namnest"/"hrrr">
def writeToCmd(stringToWrite):
    if path.exists(path.join(basePath, "plotcmds.txt")):
        currentCmdFile = open(path.join(basePath, "plotcmds.txt"), "r")
        currentStr = open(path.join(basePath, "plotcmds.txt"), "r").read()
        currentCmdFile.close()
    else:
        currentStr = ""
    if stringToWrite not in currentStr:
        with open(path.join(basePath, "plotcmds.txt"), "a") as cmdw:
            cmdw.write(stringToWrite)
            cmdw.close()

def writeToStatus(stringToWrite):
    print(stringToWrite)
    stringToWrite = stringToWrite+"\n"
    if path.exists(path.join(basePath, "status.txt")):
        currentStatusFile = open(path.join(basePath, "status.txt"), "r")
        currentStr = open(path.join(basePath, "status.txt"), "r").read()
        currentStatusFile.close()
    else:
        currentStr = ""
    if stringToWrite not in currentStr:
        with open(path.join(basePath, "status.txt"), "a") as statw:
            statw.write(stringToWrite)
            statw.close()

if __name__ == "__main__":
    basePath = path.dirname(path.abspath(__file__))
    reqVariableAddons = [
        "&lev_2_m_above_ground=on&var_TMP=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2F", #2mT
        "&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2F", # 10m Wnd
        "&lev_surface=on&var_HGT=on&var_PRES=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2F", # MSLP
        [0,1,2] #composite of 0, 1, 2
    ]
    toDate = dt(dt.utcnow().year, dt.utcnow().month, dt.utcnow().day, 0, 0, 0, 0)
    yesterDate = toDate - timedelta(days=1)
    modelName = sys.argv[1]
    writeToStatus("Beginning model fetch for model "+modelName)
    if modelName == "gfs":
        productTypeBase = 300
        initTimes = list(range(0, 19, 6))
        fHours = list(range(0, 120, 1)) + list(range(120, 385, 3))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t<MODEL_INIT_TIME>z.pgrb2.0p25.f<FHOUR_LONG><REQUESTED_VARIABLE>gfs.<MODEL_INIT_DATE>%2F<MODEL_INIT_TIME>%2Fatmos"
    elif modelName == "nam":
        productTypeBase = 500
        initTimes = list(range(0, 19, 6))
        fHours = list(range(0, 36, 1)) + list(range(36, 85, 3))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam.pl?file=nam.t<MODEL_INIT_TIME>z.awphys<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
    elif modelName == "namnest":
        productTypeBase = 600
        initTimes = list(range(0, 19, 6))
        fHours = list(range(0, 61, 1))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_nam_conusnest.pl?file=nam.t<MODEL_INIT_TIME>z.conusnest.hiresf<FHOUR_SHORT>.tm00.grib2<REQUESTED_VARIABLE>nam.<MODEL_INIT_DATE>"
    elif modelName == "hrrr":
        productTypeBase = 800
        initTimes = list(range(0, 24, 1))
        fHoursLongRun = list(range(0, 49, 1))
        fHoursShortRun = list(range(0, 19, 1))
        templateString = "https://nomads.ncep.noaa.gov/cgi-bin/filter_hrrr_2d.pl?file=hrrr.t<MODEL_INIT_TIME>z.wrfsfcf<FHOUR_SHORT>.grib2<REQUESTED_VARIABLE>hrrr.<MODEL_INIT_DATE>%2Fconus"
    else:
        raise Exception("<model> must be 'gfs', 'nam', 'namnest', or 'hrrr'")
    if path.exists(path.join(basePath, "firstPlotDT.txt")):
        readFirstPlotFile = open(path.join(basePath, "firstPlotDT.txt"), "r")
        firstPlotTime = dt.strptime(readFirstPlotFile.read(), "%Y%m%d%H%M")
        readFirstPlotFile.close()
    else:
        firstPlotTime = dt.utcnow()
        writeFirstPlotFile = open(path.join(basePath, "firstPlotDT.txt"), "w")
        writeFirstPlotFile.write(firstPlotTime.strftime("%Y%m%d%H%M"))
        writeFirstPlotFile.close()
    yesterRuns = [(yesterDate + timedelta(hours=target)).strftime("%Y%m%d%H%M") for target in initTimes]
    toRuns = [(toDate + timedelta(hours=target)).strftime("%Y%m%d%H%M") for target in initTimes if (toDate + timedelta(hours=target)) < dt.utcnow()]
    recentRuns = yesterRuns+toRuns
    for recentRun in recentRuns.copy():
        if int(recentRun) < int(firstPlotTime.strftime("%Y%m%d%H%M")):
            recentRuns.remove(recentRun)
    if modelName == "hrrr":
        hoursToRequest = dict()
        for recentRun in recentRuns:
            if int(recentRun[-4:-2]) in list(range(0, 19, 6)):
                hoursToRequest[recentRun] = fHoursLongRun
            else:
                hoursToRequest[recentRun] = fHoursShortRun
    else:
        hoursToRequest = {recentRun:fHours for recentRun in recentRuns}
    productsToRequest = [hoursToRequest for _ in range(0, len(reqVariableAddons))]
    for prodAddon in range(0, len(productsToRequest)):
        if prodAddon == 0:
            fieldToGen = "t2m"
        elif prodAddon == 1:
            fieldToGen = "sfcwind"
        elif prodAddon == 2:
            fieldToGen = "sp"
        elif prodAddon == 3:
            fieldToGen = "sfccomposite"
        productIDToCheck = prodAddon + productTypeBase
        metadataPath = path.join(basePath, "output/metadata/products/"+str(productIDToCheck)+"/")
        if path.exists(metadataPath):
            for run in recentRuns:
                runFile = path.join(metadataPath, run+".json")
                if path.exists(runFile):
                    runJsonData = dict()
                    with open(runFile) as jsonRead:
                        runJsonData = json.load(jsonRead)
                    if runJsonData["availableFrameCount"] == runJsonData["totalFrameCount"]:
                        productToTrim = productsToRequest[prodAddon].copy()
                        productToTrim.pop(run)
                        productsToRequest[prodAddon] = productToTrim
                    else:
                        frmsToDelete = [frame["fhour"] for frame in runJsonData["productFrames"]]
                        oldArr = productsToRequest[prodAddon][run]
                        newArr = [oldRun for oldRun in oldArr if oldRun not in frmsToDelete]
                        productsToRequest[prodAddon][run] = newArr
        for initRun in productsToRequest[prodAddon]:
            requestedHoursForRun = productsToRequest[prodAddon][initRun]
            variableAddonsToReq = list()
            if type(reqVariableAddons[prodAddon]) == str:
                variableAddonsToReq.append(reqVariableAddons[prodAddon])
            else:
                [variableAddonsToReq.append(reqVariableAddons[addonFromComposite]) for addonFromComposite in reqVariableAddons[prodAddon]]
            for reqVariableAddon in variableAddonsToReq:
                initDate = initRun[:-4]
                initTime = initRun[-4:-2]
                lastSuccessfulfHour = -1
                for requestedForecastHourI in requestedHoursForRun:
                    print(lastSuccessfulfHour)
                    print(requestedForecastHourI)
                    if (requestedForecastHourI - 2) > lastSuccessfulfHour:
                        print(str(requestedForecastHourI)+" is more than two fails after "+str(lastSuccessfulfHour)+"--breaking!")
                        break
                    requestedForecastHour = str(f'{requestedForecastHourI:02}')
                    requestedForecastHourLong = str(f'{requestedForecastHourI:03}')
                    shouldSkipDownload = False
                    if type(reqVariableAddons[prodAddon]) == list:
                        shouldSkipDownload = True
                        for compositeAddon in reqVariableAddons[prodAddon]:
                            if len(productsToRequest) > compositeAddon:
                                if initRun in productsToRequest[compositeAddon].keys():
                                    if requestedForecastHourI not in productsToRequest[compositeAddon][initRun]:
                                        shouldSkipDownload = False
                    if shouldSkipDownload:
                            lastSuccessfulfHour = requestedForecastHourI
                            writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun+" "+requestedForecastHour+" "+fieldToGen+"\n")
                    else:
                        urlToFetch = templateString.replace("<REQUESTED_VARIABLE>", reqVariableAddon).replace("<MODEL_INIT_TIME>", initTime).replace("<MODEL_INIT_DATE>", initDate).replace("<FHOUR_LONG>", requestedForecastHourLong).replace("<FHOUR_SHORT>", requestedForecastHour)
                        if "TMP" in reqVariableAddon:
                            saveFileName = "t2m.grib2"
                        elif "UGRD" in reqVariableAddon:
                            saveFileName = "sfcwind.grib2"
                        elif "PRES" in reqVariableAddon:
                            saveFileName = "sp.grib2"
                        modelDataPath = path.join(basePath, "modelData/"+modelName+"/"+initDate+"/"+initTime+"/"+str(requestedForecastHourI)+"/"+saveFileName)
                        Path(path.dirname(modelDataPath)).mkdir(parents=True, exist_ok=True)
                        writeToStatus("Downloading "+initRun+"Z "+modelName+" f"+requestedForecastHour+" "+fieldToGen)
                        modelData = requests.get(urlToFetch)
                        if "GRIB" in modelData.text:
                            print("download succeeded!")
                            lastSuccessfulfHour = requestedForecastHourI
                            with open(modelDataPath, "wb") as f:
                                f.write(modelData.content)
                            with open(path.join(basePath, "plotcmds.txt"), "a") as cmd:
                                writeToCmd(sys.executable+" "+path.join(basePath, "modelPlot.py")+" "+modelName+" "+initRun+" "+requestedForecastHour+" "+fieldToGen+"\n")
