#!/usr/bin/env python3
# Script for interfacing with NOAA Operational Model Archive and Distribution System (NOMADS)
# Created 17 October 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from datetime import datetime as dt, timedelta
from os import path, listdir
import json
import requests

# modelFetch.py <"gfs"/"nam"/"namnest"/"hrrr">
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
    yesterRuns = [(yesterDate + timedelta(hours=target)).strftime("%Y%m%d%H%M") for target in initTimes]
    toRuns = [(toDate + timedelta(hours=target)).strftime("%Y%m%d%H%M") for target in initTimes if (toDate + timedelta(hours=target)) < dt.utcnow()]
    recentRuns = yesterRuns+toRuns
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
                        productsToRequest[prodAddon].pop(run)
                    else:
                        frmsToDelete = [frame["fhour"] for frame in runJsonData["productFrames"]]
                        oldArr = productsToRequest[productIDToCheck][run]
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
                for requestedForecastHourI in requestedHoursForRun:
                    shouldSkipDownload = False
                    if type(reqVariableAddons[prodAddon]) == list:
                        shouldSkipDownload = True
                        for compositeAddon in reqVariableAddons[prodAddon]:
                            if requestedForecastHourI not in productsToRequest[compositeAddon][initRun]:
                                shouldSkipDownload = False
                    if shouldSkipDownload:
                        logfile = open("debug.txt", "a")
                        logfile.write("Skipped DL...\n")
                        logfile.close()
                    else:
                        requestedForecastHour = str(f'{requestedForecastHourI:02}')
                        requestedForecastHourLong = str(f'{requestedForecastHourI:03}')
                        urlToFetch = templateString.replace("<REQUESTED_VARIABLE>", reqVariableAddon).replace("<MODEL_INIT_TIME>", initTime).replace("<MODEL_INIT_DATE>", initDate).replace("<FHOUR_LONG>", requestedForecastHourLong).replace("<FHOUR_SHORT>", requestedForecastHour)
                    if "GRIB" in requests.get(urlToFetch).text:
                        logfile = open("debug.txt", "a")
                        logfile.write(".")
                        logfile.close()
                    else:
                        logfile = open("debug.txt", "a")
                        logfile.write("\n\n")
                        logfile.write(modelName+"\n")
                        logfile.write(initDate+"\n")
                        logfile.write(initTime+"\n")
                        logfile.write(requestedForecastHourLong+"\n")
                        logfile.write(str(prodAddon))
                        logfile.write("\n\n")
                        logfile.close()
                    