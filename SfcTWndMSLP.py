#!/usr/bin/env python3
# Surface temp, wind, mslp model postprocessing for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path, remove, listdir
from pathlib import Path
import xarray as xr
from metpy import constants
from metpy.plots import USCOUNTIES
import metpy
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime as dt, timedelta
from pandas import Timestamp as pdtimestamp
from matplotlib import image as mpimage
from scipy import ndimage
import json
import warnings

# SfcTWndMSLP.py <model> <initialization> <fhour>
modelName = sys.argv[1]
initDateTime = dt.strptime(sys.argv[2], "%Y%m%d%H%M")
fhour = int(sys.argv[3])
basePath = path.dirname(path.abspath(__file__))
axExtent = [-130, -60, 20, 50]
if modelName == "gfs":
    productTypeBase = 300
    productFrameCount = 209
elif modelName == "nam":
    productTypeBase = 500
    productFrameCount = 53
elif modelName == "namnest":
    productTypeBase = 600
    productFrameCount = 61
elif modelName == "hrrr":
    productTypeBase = 800
    if initDateTime.hour in [0, 6, 12, 18]:
        productFrameCount = 49
    else:
        productFrameCount = 19
else:
    raise Exception("<model> must be 'gfs', 'nam', 'namnest', or 'hrrr'")

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


def writeJson(productID, gisInfo):
    if productID == 300:
        productDesc = "GFS Surface Temperature"
        dirname = "sfcT"
    elif productID == 301:
        productDesc = "GFS Surface Winds"
        dirname = "sfcWnd"
    elif productID == 302:
        productDesc = "GFS Surface MSLP"
        dirname = "sfcMSLP"
    elif productID == 303:
        productDesc = "GFS Surface Temperature, Winds, MSLP"
        dirname = "sfcTWndMSLP"
    elif productID == 500:
        productDesc = "NAM Surface Temperature"
        dirname = "sfcT"
    elif productID == 501:
        productDesc = "NAM Surface Winds"
        dirname = "sfcWnd"
    elif productID == 502:
        productDesc = "NAM Surface MSLP"
        dirname = "sfcMSLP"
    elif productID == 503:
        productDesc = "NAM Surface Temperature, Winds, MSLP"
        dirname = "sfcTWndMSLP"
    elif productID == 600:
        productDesc = "NAM NEST Surface Temperature"
        dirname = "sfcT"
    elif productID == 601:
        productDesc = "NAM NEST Surface Winds"
        dirname = "sfcWnd"
    elif productID == 602:
        productDesc = "NAM NEST Surface MSLP"
        dirname = "sfcMSLP"
    elif productID == 603:
        productDesc = "NAM NEST Surface Temperature, Winds, MSLP"
        dirname = "sfcTWndMSLP"
    elif productID == 800:
        productDesc = "HRRR Surface Temperature"
        dirname = "sfcT"
    elif productID == 801:
        productDesc = "HRRR Surface Winds"
        dirname = "sfcWnd"
    elif productID == 802:
        productDesc = "HRRR Surface MSLP"
        dirname = "sfcMSLP"
    elif productID == 803:
        productDesc = "HRRR Surface Temperature, Winds, MSLP"
        dirname = "sfcTWndMSLP"
    if gisInfo == ["0,0", "0,0"]:
        isGIS = False
        productPath = "products/"
    else:
        isGIS = True
        productPath = "gisproducts/"
    productPath = productPath+modelName+"/"+dirname+"/"
    pathExtension = initDateTime.strftime("%Y/%m/%d/%H%M/")
    publishTime = dt.utcnow()
    productDict = {
        "productID" : productID,
        "productDescription" : productDesc,
        "productPath" : productPath,
        "productReloadTime" : 300,
        "lastReloadTime" : publishTime.strftime("%Y%m%d%H%M"),
        "isForecast" : True,
        "isGIS" : isGIS
    }
    productDictJsonPath = path.join(basePath, "output/metadata/"+str(productID)+".json")
    Path(path.dirname(productDictJsonPath)).mkdir(parents=True, exist_ok=True)
    if path.exists(productDictJsonPath):
        remove(productDictJsonPath)
    with open(productDictJsonPath, "w") as jsonWrite:
        json.dump(productDict, jsonWrite, indent=4)
    productRunDictPath = path.join(basePath, "output/metadata/products/"+str(productID)+"/"+initDateTime.strftime("%Y%m%d%H%M")+".json")
    Path(path.dirname(productRunDictPath)).mkdir(parents=True, exist_ok=True)
    if path.exists(productRunDictPath):
        remove(productRunDictPath)
    framesArray = list()
    if len(listdir(path.join(basePath, "output/"+productPath+pathExtension))) > 0:
        frameNames = listdir(path.join(basePath, "output/"+productPath+pathExtension))
        frameHours = [int(framename.replace("f", "").replace(".png", "")) for framename in frameNames]
        for frameHr in frameHours:
            fvalidTime = initDateTime + timedelta(hours=frameHr)
            frmDict = {
                "fhour" : frameHr,
                "filename" : "f"+str(frameHr)+".png",
                "gisInfo" : gisInfo,
                "valid" : int(fvalidTime.strftime("%Y%m%d%H%M"))
            }
            framesArray.append(frmDict)
    productRunDict = {
        "publishTime" : publishTime.strftime("%Y%m%d%H%M"),
        "pathExtension" : pathExtension,
        "runName" : initDateTime.strftime("%d %b %Y %HZ"),
        "availableFrameCount" : len(framesArray),
        "totalFrameCount" : productFrameCount,
        "productFrames" : sorted(framesArray, key=lambda dict: dict["fhour"]) 
    }
    with open(productRunDictPath, "w") as jsonWrite:
        json.dump(productRunDict, jsonWrite, indent=4)
    

def set_size(w,h, ax=None):
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

def staticSFCTempWindMSLPPlot(tempArr, wind, mslp):
    temp = tempArr[0]
    contourHandle = tempArr[1]
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    set_size(1920*px, 1080*px, ax=ax)
    tempImg = mpimage.imread(temp)
    ax.imshow(tempImg, extent=axExtent, transform=ccrs.PlateCarree(), zorder=1)
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5, zorder=3)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5, zorder=3)
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    fig.colorbar(contourHandle, cax=cbax, orientation="horizontal")
    cbax.set_xlabel("Temperature (°F)")
    wndImg = mpimage.imread(wind)
    ax.imshow(wndImg, extent=axExtent, transform=ccrs.PlateCarree(), zorder=5)
    mslpImg = mpimage.imread(mslp)
    ax.imshow(mslpImg, extent=axExtent, transform=ccrs.PlateCarree(), zorder=6)
    validTime = initDateTime + timedelta(hours=fhour)
    tax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+.01,0.045,(ax.get_position().width/3),.05])
    tax.text(0.5, 0.5, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n2m Temp, 10m Winds, MSLP\nf"+str(fhour)+" Valid "+validTime.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
    tax.set_xlabel("Python HDWX -- Send bugs to stgardner4@tamu.edu")
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    lax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+tax.get_position().width+.01,0,(ax.get_position().width/3),.1])
    lax.set_aspect(2821/11071)
    lax.axis("off")
    plt.setp(lax.spines.values(), visible=False)
    atmoLogo = mpimage.imread("assets/atmoLogo.png")
    lax.imshow(atmoLogo)
    fig.set_facecolor("white")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output/products/"+modelName+"/sfcTWndMSLP/"+runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig.savefig(staticSavePath+"/f"+str(fhour)+".png", bbox_inches="tight")
    plt.close(fig)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 3
    writeJson(productId, gisInfo)

def tempPlot(pathToRead):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
        # modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib", filter_by_keys={"typeOfLevel": "heightAboveGround", "topLevel": 2})
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcT/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    gisSavePath = path.join(gisSavePath, "f"+str(fhour)+".png")
    tempData = modelDataArray["t"]
    tempData = tempData.metpy.quantify()
    tempData = tempData.metpy.convert_units("degF")
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs":
        shouldTransformFirst = False
    else:
        shouldTransformFirst = True
    contourmap = ax.contourf(tempData.longitude, tempData.latitude, tempData, levels=np.arange(-20, 120, 5), cmap="nipy_spectral", vmin=-20, vmax=120, transform=ccrs.PlateCarree(), transform_first=shouldTransformFirst)
    set_size(1920*px, 1080*px, ax=ax)
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    plt.close(fig)
    gisInfo = ["20,-130", "50,-60"]
    writeJson(productTypeBase, gisInfo)
    return [gisSavePath, contourmap]

def windPlot(pathToRead):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
        # modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib", filter_by_keys={"typeOfLevel": "heightAboveGround", "topLevel": 10})
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcWnd/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    gisSavePath = path.join(gisSavePath, "f"+str(fhour)+".png")
    uwnd = modelDataArray["u10"]
    uwnd = uwnd.metpy.quantify()
    uwnd = uwnd.metpy.convert_units("kt")
    vwnd = modelDataArray["v10"]
    vwnd = vwnd.metpy.quantify()
    vwnd = vwnd.metpy.convert_units("kt")
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs":
        spatialLimit = slice(None, None, 5)
        dataLimit = (slice(None, None, 5), slice(None, None, 5))
    elif modelName == "nam":
        spatialLimit = (slice(None, None, 10), slice(None, None, 10))
        dataLimit = (slice(None, None, 10), slice(None, None, 10))
    elif modelName == "namnest" or modelName == "hrrr":
        spatialLimit = (slice(None, None, 40), slice(None, None, 40))
        dataLimit = (slice(None, None, 40), slice(None, None, 40))
    windbarbs = ax.barbs(uwnd.longitude.data[spatialLimit], uwnd.latitude.data[spatialLimit], uwnd.data[dataLimit], vwnd.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5)
    set_size(1920*px, 1080*px, ax=ax)
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    plt.close(fig)
    gisInfo = ["20,-130", "50,-60"]
    productId = productTypeBase + 1
    writeJson(productId, gisInfo)
    return gisSavePath

def mslpPlot(pathToRead):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
        # modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib", filter_by_keys={"typeOfLevel": "surface"})
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    gisSavePath = path.join(gisSavePath, "f"+str(fhour)+".png")
    barometricPressData = modelDataArray["sp"]
    barometricPressData = barometricPressData.metpy.quantify()
    orogData = modelDataArray["orog"]
    orogData = orogData.metpy.quantify()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        # tempData = xr.open_dataset(pathToRead, engine="cfgrib", filter_by_keys={"typeOfLevel": "heightAboveGround", "topLevel": 2})
    tempData = modelDataArray["t"]
    tempData = tempData.metpy.quantify()
    # I tried using mpcalc altimeter->mslp function here, but it ended up doing nothing and I don't feel like figuring out why
    # Therefore I implemented the same equation manually...
    mslpData = barometricPressData * np.exp(orogData*constants.earth_gravity/(constants.dry_air_gas_constant*tempData))
    mslpData = mslpData.metpy.quantify()
    mslpData = mslpData.metpy.convert_units("hPa")
    # Unidata says smoothing MSLP "a little" is... well they didn't comment on why, they just did it, and it makes the rocky mtns less noisy...
    # https://unidata.github.io/python-gallery/examples/MSLP_temp_winds.html
    if modelName == "namnest" or modelName == "hrrr":
        mslpData = ndimage.gaussian_filter(mslpData.data, 5)
    else:
        mslpData = ndimage.gaussian_filter(mslpData.data, 3)
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs":
        shouldTransformFirst = False
    else:
        shouldTransformFirst = True
    contourmap = ax.contour(barometricPressData.longitude, barometricPressData.latitude, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree(), transform_first=shouldTransformFirst)
    contourLabels = ax.clabel(contourmap, levels=np.arange(800, 1200, 2), inline=True, fontsize=15)
    [label.set_rotation(0) for label in contourLabels]    
    set_size(1920*px, 1080*px, ax=ax)
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(gisSavePath, transparent=True, bbox_inches=extent)
    plt.close(fig)
    gisInfo = ["20,-130", "50,-60"]
    productId = productTypeBase + 2
    writeJson(productId, gisInfo)
    return gisSavePath

if __name__ == "__main__":
    writeToStatus(str("Plotting init "+str(initDateTime.hour)+"Z f"+str(fhour)+" "+modelName+" sfcTWndMSLP"))
    inputPath = path.join(basePath, "modelData/"+modelName+"/"+dt.strftime(initDateTime, "%Y%m%d")+"/"+dt.strftime(initDateTime, "%H")+"/"+str(fhour))
    pathToRead = path.join(inputPath, "sfcTWndMSLP.grib2")
    if path.exists(pathToRead):
        staticSFCTempWindMSLPPlot(tempPlot(pathToRead), windPlot(pathToRead), mslpPlot(pathToRead))
    [remove(path.join(inputPath, idxFile)) for idxFile in listdir(inputPath) if idxFile.endswith("idx")]