#!/usr/bin/env python3
# Surface temp, wind, mslp model postprocessing for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path, listdir
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

# modelPlot.py <model> <initialization> <fhour> <field to plot>
modelName = sys.argv[1]
initDateTime = dt.strptime(sys.argv[2], "%Y%m%d%H%M")
fhour = int(sys.argv[3])
fieldToPlot = sys.argv[4]
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
        "productReloadTime" : 60,
        "lastReloadTime" : publishTime.strftime("%Y%m%d%H%M"),
        "isForecast" : True,
        "isGIS" : isGIS,
        "fileExtension" : "png"
    }
    productDictJsonPath = path.join(basePath, "output/metadata/"+str(productID)+".json")
    Path(path.dirname(productDictJsonPath)).mkdir(parents=True, exist_ok=True)
    with open(productDictJsonPath, "w") as jsonWrite:
        json.dump(productDict, jsonWrite, indent=4)
    productRunDictPath = path.join(basePath, "output/metadata/products/"+str(productID)+"/"+initDateTime.strftime("%Y%m%d%H%M")+".json")
    Path(path.dirname(productRunDictPath)).mkdir(parents=True, exist_ok=True)
    if path.exists(productRunDictPath):
        with open(productRunDictPath, "r") as jsonRead:
            oldData = json.load(jsonRead)
        framesArray = oldData["productFrames"]
    else:
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
            if frmDict not in framesArray:
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
    productTypeID = int(str(productTypeBase)[0])
    productTypeDictPath = path.join(basePath, "output/metadata/productTypes/"+str(productTypeID)+".json")
    Path(path.dirname(productTypeDictPath)).mkdir(parents=True, exist_ok=True)
    productsInType = list()
    if path.exists(productTypeDictPath):
        with open(productTypeDictPath, "r") as jsonRead:
            oldProductTypeDict = json.load(jsonRead)
        for productInOldDict in oldProductTypeDict["products"]:
            if productInOldDict["productID"] != productID:
                productsInType.append(productInOldDict)
    productsInType.append(productDict)
    productTypeDict = {
        "productTypeID" : productTypeID,
        "productTypeDescription" : modelName.upper(),
        "products" : sorted(productsInType, key=lambda dict: dict["productID"])
    }
    with open(productTypeDictPath, "w") as jsonWrite:
        json.dump(productTypeDict, jsonWrite, indent=4)
    

def set_size(w,h, ax=None):
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

def staticSFCTempWindMSLPPlot():
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-131, -61, 21, 53], crs=ccrs.PlateCarree())
    set_size(1920*px, 1080*px, ax=ax)
    contourmap = tempPlot(False, ax=ax)
    windPlot(False, ax=ax)
    mslpPlot(False, ax=ax)
    ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal")
    cbax.set_xlabel("Temperature (°F)")
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

def tempPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "t2m.grib2")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcT/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    tempData = modelDataArray["t2m"]
    tempData = tempData.metpy.quantify()
    tempData = tempData.metpy.convert_units("degF")
    validTime = pdtimestamp(np.datetime64(tempData.valid_time.data)).to_pydatetime()
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs":
        lonsToPlot = np.tile(np.array([tempData.longitude.data]), (tempData.data.shape[0], 1))
        latsToPlot = np.tile(tempData.latitude.data, (tempData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = tempData.longitude
        latsToPlot = tempData.latitude
    contourmap = ax.contourf(lonsToPlot, latsToPlot, tempData, levels=np.arange(-20, 120, 5), cmap="nipy_spectral", vmin=-20, vmax=120, transform=ccrs.PlateCarree(), transform_first=True)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        writeJson(productTypeBase, gisInfo)
    return contourmap

def windPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sfcwind.grib2")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcWnd/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    uwind = modelDataArray["u10"]
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("kt")
    vwind = modelDataArray["v10"]
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("kt")
    validTime = pdtimestamp(np.datetime64(modelDataArray.valid_time.data)).to_pydatetime()
    if standaloneFig:
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
    windbarbs = ax.barbs(uwind.longitude.data[spatialLimit], uwind.latitude.data[spatialLimit], uwind.data[dataLimit], vwind.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 1
        writeJson(productId, gisInfo)
    return windbarbs

def mslpPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sp.grib2")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    barometricPressData = modelDataArray["sp"]
    barometricPressData = barometricPressData.metpy.quantify()
    orogData = modelDataArray["orog"]
    orogData = orogData.metpy.quantify()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        tempData = xr.open_dataset(path.join(inputPath, "t2m.grib2"), engine="cfgrib")
    tempData = tempData["t2m"]
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
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs":
        lonsToPlot = np.tile(np.array([barometricPressData.longitude.data]), (barometricPressData.data.shape[0], 1))
        latsToPlot = np.tile(barometricPressData.latitude.data, (barometricPressData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = barometricPressData.longitude
        latsToPlot = barometricPressData.latitude
    contourmap = ax.contour(lonsToPlot, latsToPlot, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree(), transform_first=True)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    contourLabels = ax.clabel(contourmap, levels=np.arange(800, 1200, 2), inline=True, fontsize=15)
    [label.set_rotation(0) for label in contourLabels]
    if standaloneFig:
        validTime = pdtimestamp(np.datetime64(modelDataArray.valid_time.data)).to_pydatetime()
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 2
        writeJson(productId, gisInfo)
    return contourmap

if __name__ == "__main__":
    writeToStatus(str("Plotting init "+str(initDateTime.hour)+"Z f"+str(fhour)+" "+modelName+" "+fieldToPlot))
    inputPath = path.join(basePath, "modelData/"+modelName+"/"+dt.strftime(initDateTime, "%Y%m%d")+"/"+dt.strftime(initDateTime, "%H")+"/"+str(fhour))
    sfcTempPath = path.join(inputPath, "t2m.grib2")
    sfcWindsPath = path.join(inputPath, "sfcwind.grib2")
    sfcPressPath = path.join(inputPath, "sp.grib2")
    if path.exists(inputPath):
        if fieldToPlot == "t2m" and path.exists(sfcTempPath):
            tempPlot(True)
        if fieldToPlot == "sfcwind" and path.exists(sfcWindsPath):
            windPlot(True)
        if fieldToPlot == "sp" and path.exists(sfcPressPath) and path.exists(sfcTempPath):
            mslpPlot(True)
        if fieldToPlot == "sfccomposite" and path.exists(sfcTempPath) and path.exists(sfcWindsPath) and path.exists(sfcPressPath):
            staticSFCTempWindMSLPPlot()
