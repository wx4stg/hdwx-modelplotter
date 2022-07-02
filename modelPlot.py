#!/usr/bin/env python3
# Model postprocessing for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path, listdir, remove, chmod
from pathlib import Path
import xarray as xr
from metpy import constants
from metpy.plots import USCOUNTIES
from metpy import calc as mpcalc
import metpy
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime as dt, timedelta
from pandas import Timestamp as pdtimestamp
from matplotlib import image as mpimage
from matplotlib import colors
from scipy import ndimage
import json
import warnings
from atomicwrites import atomic_write

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
elif modelName == "ecmwf-hres":
    productTypeBase = 1000
    if initDateTime.hour in [0, 12]:
        productFrameCount = 61
    else:
        productFrameCount = 31
else:
    raise Exception("<model> must be 'gfs', 'nam', 'namnest', 'hrrr', or 'ecmwf-hres'")

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
    elif productID == 316:
        productDesc = "GFS 500 hPa Winds"
        dirname = "500wind"
    elif productID == 321:
        productDesc = "GFS 250 hPa Winds"
        dirname = "250wind"
    elif productID == 325:
        productDesc = "GFS 850 hPa Winds"
        dirname = "850wind"
    elif productID == 390:
        productDesc = "GFS Surface Wind Divergence"
        dirname = "divergence"
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
    elif productID == 516:
        productDesc = "NAM 500 hPa Winds"
        dirname = "500wind"
    elif productID == 521:
        productDesc = "NAM 250 hPa Winds"
        dirname = "250wind"
    elif productID == 525:
        productDesc = "NAM 850 hPa Winds"
        dirname = "850wind"
    elif productID == 590:
        productDesc = "NAM Surface Wind Divergence"
        dirname = "divergence"
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
    elif productID == 616:
        productDesc = "NAM NEST 500 hPa Winds"
        dirname = "500wind"
    elif productID == 621:
        productDesc = "NAM NEST 250 hPa Winds"
        dirname = "250wind"
    elif productID == 625:
        productDesc = "NAM NEST 850 hPa Winds"
        dirname = "850wind"
    elif productID == 690:
        productDesc = "NAM NEST Surface Wind Divergence"
        dirname = "divergence"
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
    elif productID == 816:
        productDesc = "HRRR 500 hPa Winds"
        dirname = "500wind"
    elif productID == 821:
        productDesc = "HRRR 250 hPa Winds"
        dirname = "250wind"
    elif productID == 825:
        productDesc = "HRRR 850 hPa Winds"
        dirname = "850wind"
    elif productID == 890:
        productDesc = "HRRR Surface Wind Divergence"
        dirname = "divergence"
    elif productID == 1000:
        productDesc = "ECMWF-HRES Surface Temperature"
        dirname = "sfcT"
    elif productID == 1001:
        productDesc = "ECMWF-HRES Surface Winds"
        dirname = "sfcWnd"
    elif productID == 1002:
        productDesc = "ECMWF-HRES Surface MSLP"
        dirname = "sfcMSLP"
    elif productID == 1003:
        productDesc = "ECMWF-HRES Surface Temperature, Winds, MSLP"
        dirname = "sfcTWndMSLP"
    elif productID == 1016:
        productDesc = "ECMWF-HRES 500 hPa Winds"
        dirname = "500wind"
    elif productID == 1021:
        productDesc = "ECMWF-HRES 250 hPa Winds"
        dirname = "250wind"
    elif productID == 1025:
        productDesc = "ECMWF-HRES 850 hPa Winds"
        dirname = "850wind"
    elif productID == 1090:
        productDesc = "ECMWF Surface Wind Divergence"
        dirname = "divergence"
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
        "fileExtension" : "png",
        "displayFrames" : 0
    }
    productDictJsonPath = path.join(basePath, "output/metadata/"+str(productID)+".json")
    Path(path.dirname(productDictJsonPath)).mkdir(parents=True, exist_ok=True)
    with atomic_write(productDictJsonPath, overwrite=True) as jsonWrite:
        json.dump(productDict, jsonWrite, indent=4)
    chmod(productDictJsonPath, 0o644)
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
        frameHours = [int(framename.replace("f", "").replace(".png", "")) for framename in frameNames if ".png" in framename]
        for frameHr in frameHours:
            fvalidTime = initDateTime + timedelta(hours=frameHr)
            frmDict = {
                "fhour" : frameHr,
                "filename" : "f"+str(frameHr)+".png",
                "gisInfo" : gisInfo,
                "valid" : fvalidTime.strftime("%Y%m%d%H%M")
            }
            if frmDict not in framesArray:
                framesArray.append(frmDict)
    productRunDict = {
        "publishTime" : int(publishTime.strftime("%Y%m%d%H%M")),
        "pathExtension" : pathExtension,
        "runName" : initDateTime.strftime("%d %b %Y %HZ"),
        "availableFrameCount" : len(framesArray),
        "totalFrameCount" : productFrameCount,
        "productFrames" : sorted(framesArray, key=lambda dict: dict["fhour"]) 
    }
    with atomic_write(productRunDictPath, overwrite=True) as jsonWrite:
        json.dump(productRunDict, jsonWrite, indent=4)
    chmod(productRunDictPath, 0o644)
    productTypeID = productTypeBase // 100
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
    with atomic_write(productTypeDictPath, overwrite=True) as jsonWrite:
        json.dump(productTypeDict, jsonWrite, indent=4)
    chmod(productTypeDictPath, 0o644)

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
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    contourmap = sfcTempPlot(False, ax=ax)
    sfcWindPlot(False, ax=ax)
    mslpPlot(False, ax=ax)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal", extend="both").set_ticks(np.sort(np.append(np.arange(-40, 120, 10), 32)))
    targetIdx = list(list(plt.xticks())[0]).index(32)
    plt.xticks()[-1][targetIdx].set_color("red")
    cbax.set_xlabel("Temperature (°F)")
    validTime = initDateTime + timedelta(hours=fhour)
    fig.set_size_inches(1920*px, 1080*px)
    tax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+.01,0.045,(ax.get_position().width/3),.05])
    tax.text(0.5, 0.5, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n2m Temp, 10m Winds, MSLP\nf"+str(fhour)+" Valid "+validTime.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
    if "ecmwf" in modelName:
        tax.set_xlabel("\nPython HDWX -- Send bugs to stgardner4@tamu.edu\nCopyright © "+initDateTime.strftime("%Y")+" European Centre for Medium-Range Weather Forecasts (ECMWF)\nhttps://www.ecmwf.int/")
    else:
        tax.set_xlabel("\nPython HDWX -- Send bugs to stgardner4@tamu.edu\nData provided by NOAA/NCEP")
    plt.setp(tax.spines.values(), visible=False)
    tax.tick_params(left=False, labelleft=False)
    tax.tick_params(bottom=False, labelbottom=False)
    lax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+tax.get_position().width+.01,0,(ax.get_position().width/3),.1])
    lax.set_aspect(2821/11071)
    lax.axis("off")
    plt.setp(lax.spines.values(), visible=False)
    atmoLogo = mpimage.imread(path.join(basePath, "assets", "atmoLogo.png"))
    lax.imshow(atmoLogo)
    ax.set_position([ax.get_position().x0, .12, ax.get_position().width, ax.get_position().height])
    fig.set_facecolor("white")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output/products/"+modelName+"/sfcTWndMSLP/"+runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig.savefig(staticSavePath+"/f"+str(fhour)+".png", bbox_inches="tight")
    plt.close(fig)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 3
    writeJson(productId, gisInfo)

def sfcTempPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "t2m.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcT/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    tempData = modelDataArray["t2m"]
    tempData = tempData.metpy.quantify()
    tempData = tempData.metpy.convert_units("degF")
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([tempData.longitude.data]), (tempData.data.shape[0], 1))
        latsToPlot = np.tile(tempData.latitude.data, (tempData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = tempData.longitude
        latsToPlot = tempData.latitude
    frozenColorMap = plt.cm.cool_r(np.linspace(0, 1, 25))
    meltedColorMap = plt.cm.turbo(np.linspace(0.35, 1, 25))
    all_colors = np.vstack((frozenColorMap, meltedColorMap))
    temperatureColorMap = colors.LinearSegmentedColormap.from_list("temperatureColorMap", all_colors)
    temperatureNorm = colors.TwoSlopeNorm(vcenter=32, vmin=-40, vmax=130)
    levelsToContour = np.sort(np.append(np.arange(-40, 120, 10), 32))
    contourmap = ax.contourf(lonsToPlot, latsToPlot, tempData, levels=levelsToContour,  cmap=temperatureColorMap, norm=temperatureNorm, transform=ccrs.PlateCarree(), transform_first=True)
    ax.contour(lonsToPlot, latsToPlot, tempData, levels=[32], colors="red", transform=ccrs.PlateCarree(), transform_first=True, linewidths=1)
    contourLabels = ax.clabel(contourmap, levels=np.arange(-40, 120, 20), inline=True, fontsize=10, colors="black")
    [label.set_rotation(0) for label in contourLabels]
    [label.set_text(label.get_text()+"F") for label in contourLabels]
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        writeJson(productTypeBase, gisInfo)
    return contourmap

def sfcWindPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sfcwind.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output", "gisproducts", modelName, "sfcWnd", runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    uwind = modelDataArray["u10"]
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("kt")
    vwind = modelDataArray["v10"]
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("kt")
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        spatialLimit = slice(None, None, 5)
        dataLimit = (slice(None, None, 5), slice(None, None, 5))
    elif modelName == "nam":
        spatialLimit = (slice(None, None, 10), slice(None, None, 10))
        dataLimit = (slice(None, None, 10), slice(None, None, 10))
    elif modelName == "namnest" or modelName == "hrrr":
        spatialLimit = (slice(None, None, 40), slice(None, None, 40))
        dataLimit = (slice(None, None, 40), slice(None, None, 40))
    windbarbs = ax.barbs(uwind.longitude.data[spatialLimit], uwind.latitude.data[spatialLimit], uwind.data[dataLimit], vwind.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5, linewidth=0.5, zorder=2)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 1
        writeJson(productId, gisInfo)
        
        # TRACER sea-breeze convergence thing
        div = mpcalc.divergence(uwind.metpy.convert_units("m/s"), vwind.metpy.convert_units("m/s"))
        ax.scatter([-95.74010, -95.059671], [29.32937, 29.669285], s=50, c="green", marker="s", transform=ccrs.PlateCarree(), zorder=3)
        contourmap = ax.contourf(div.longitude, div.latitude, div, cmap="Reds_r", levels=np.linspace(-.004, 0, 100), transform=ccrs.PlateCarree(), transform_first=True, zorder=1)
        spatialLimit = (slice(None, None, 10), slice(None, None, 10))
        dataLimit = (slice(None, None, 10), slice(None, None, 10))
        windbarbs = ax.barbs(uwind.longitude.data[spatialLimit], uwind.latitude.data[spatialLimit], uwind.data[dataLimit], vwind.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5, linewidth=0.5, zorder=2)
        tracerAxExtent = [-99.5, -91, 26, 33.5]
        ax.set_extent(tracerAxExtent)
        # For the "static"/non-GIS/opaque image, add county/state/coastline borders
        ax.add_feature(USCOUNTIES.with_scale("5m"), edgecolor="gray", zorder=2)
        ax.add_feature(cfeat.STATES.with_scale("10m"), linewidth=0.5, zorder=3)
        ax.add_feature(cfeat.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)
        # Reduce whitespace around data axes
        ax.set_box_aspect(9/16)
        # Move the data axes to maximize the amount of space available to it
        ax.set_position([0.05, 0.11, .9, .87])
        cbax = fig.add_axes([.01,0.075,(ax.get_position().width/3),.02])
        cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal", label="Divergence (1/s)", extend="min")
        cb.set_ticks(np.arange(-.004, 0, .001))
        cb.formatter.set_powerlimits((0,0))
        tax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+.01,0.03,(ax.get_position().width/3),.05])
        validTime = initDateTime + timedelta(hours=fhour)
        tax.text(0.5, 0.5, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\nSurface wind convergence\nf"+str(fhour)+" Valid "+validTime.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
        tax.set_xlabel("Python HDWX -- Send bugs to stgardner4@tamu.edu")
        plt.setp(tax.spines.values(), visible=False)
        tax.tick_params(left=False, labelleft=False)
        tax.tick_params(bottom=False, labelbottom=False)
        lax = fig.add_axes([(.99-(ax.get_position().width/3)),0,(ax.get_position().width/3),.1])
        lax.set_aspect(2821/11071)
        lax.axis("off")
        plt.setp(lax.spines.values(), visible=False)
        atmoLogo = mpimage.imread(path.join(basePath, "assets", "atmoLogo.png"))
        lax.imshow(atmoLogo)
        ax.set_position([.005, cbax.get_position().y0+cbax.get_position().height+.005, .99, (.99-(cbax.get_position().y0+cbax.get_position().height))])
        # Set size to 1080p, resolution of the weather center monitors
        fig.set_size_inches(1920*px, 1080*px)
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        staticSavePath = path.join(basePath, "output/products/"+modelName+"/divergence/"+runPathExt)
        # Create save directory if it doesn't already exist
        Path(staticSavePath).mkdir(parents=True, exist_ok=True)
        # Write the image
        fig.savefig(staticSavePath+"/f"+str(fhour)+".png")
        writeJson(productTypeBase + 90, ["0,0", "0,0"])
        plt.close(fig)
        return windbarbs

def mslpPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sp.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    if modelName == "ecmwf-hres":
        mslpData = modelDataArray["msl"]
        mslpData = mslpData.metpy.quantify()
        mslpData = mslpData.metpy.convert_units("hPa")
    else:
        barometricPressData = modelDataArray["sp"]
        barometricPressData = barometricPressData.metpy.quantify()
        orogData = modelDataArray["orog"]
        orogData = orogData.metpy.quantify()
        [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
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
        mslpData.data = ndimage.gaussian_filter(mslpData.data, 5)
    else:
        mslpData.data = ndimage.gaussian_filter(mslpData.data, 3)
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([mslpData.longitude.data]), (mslpData.data.shape[0], 1))
        latsToPlot = np.tile(mslpData.latitude.data, (mslpData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = mslpData.longitude
        latsToPlot = mslpData.latitude
    contourmap = ax.contour(lonsToPlot, latsToPlot, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree(), transform_first=True, linewidths=0.5)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    contourLabels = ax.clabel(contourmap, levels=np.arange(800, 1040, 2), inline=True, fontsize=10)
    [label.set_rotation(0) for label in contourLabels]
    if standaloneFig:
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 2
        writeJson(productId, gisInfo)
    return contourmap

def windsAtHeightPlot(pressureLevel, standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "winds.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"wind", runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    uwind = modelDataArray["u"]
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("kt")
    vwind = modelDataArray["v"]
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("kt")
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        spatialLimit = slice(None, None, 5)
        dataLimit = (slice(None, None, 5), slice(None, None, 5))
    elif modelName == "nam":
        spatialLimit = (slice(None, None, 10), slice(None, None, 10))
        dataLimit = (slice(None, None, 10), slice(None, None, 10))
    elif modelName == "namnest" or modelName == "hrrr":
        spatialLimit = (slice(None, None, 40), slice(None, None, 40))
        dataLimit = (slice(None, None, 40), slice(None, None, 40))
    windbarbs = ax.barbs(uwind.longitude.data[spatialLimit], uwind.latitude.data[spatialLimit], uwind.data[dataLimit], vwind.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5, linewidth=0.5)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if pressureLevel == 250:
            addon = 21
        elif pressureLevel == 500:
            addon = 16
        elif pressureLevel == 850:
            addon = 25
        productId = productTypeBase + addon
        writeJson(productId, gisInfo)
    return windbarbs
if __name__ == "__main__":
    writeToStatus(str("Plotting init "+str(initDateTime.hour)+"Z f"+str(fhour)+" "+modelName+" "+fieldToPlot))
    inputPath = path.join(basePath, "modelData/"+modelName+"/"+dt.strftime(initDateTime, "%Y%m%d")+"/"+dt.strftime(initDateTime, "%H")+"/"+str(fhour))
    sfcTempPath = path.join(inputPath, "t2m.grib2")
    sfcWindsPath = path.join(inputPath, "sfcwind.grib2")
    sfcPressPath = path.join(inputPath, "sp.grib2")
    windsAtHeightPath = path.join(inputPath, "winds.grib2")
    if path.exists(inputPath):
        if fieldToPlot == "t2m" and path.exists(sfcTempPath):
            sfcTempPlot(True)
        if fieldToPlot == "sfcwind" and path.exists(sfcWindsPath):
            sfcWindPlot(True)
        if fieldToPlot == "sp" and path.exists(sfcPressPath) and path.exists(sfcTempPath):
            mslpPlot(True)
        if fieldToPlot == "sfccomposite" and path.exists(sfcTempPath) and path.exists(sfcWindsPath) and path.exists(sfcPressPath):
            staticSFCTempWindMSLPPlot()
        if fieldToPlot == "winds" and path.exists(windsAtHeightPath):
            [windsAtHeightPlot(pressSfc, True) for pressSfc in [250, 500, 850]]
