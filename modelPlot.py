#!/usr/bin/env python3
# Model postprocessing for python-based HDWX
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
from matplotlib import colors as pltcolors
from scipy import ndimage
import json
import warnings
from atomicwrites import atomic_write

# modelPlot.py <model> <initialization> <fhour> <field to plot>
modelName = sys.argv[1]
initDateTime = dt.strptime(sys.argv[2], "%Y%m%d%H%M")
fhour = int(sys.argv[3])
fieldToPlot = sys.argv[4]
axExtent = [-130, -60, 20, 50]
basePath = path.dirname(path.abspath(__file__))

hasHelpers = False
if path.exists(path.join(basePath, "HDWX_helpers.py")):
    import HDWX_helpers
    hasHelpers = True

if modelName == "gfs":
    productTypeBase = 300
elif modelName == "nam":
    productTypeBase = 500
elif modelName == "namnest":
    productTypeBase = 600
elif modelName == "hrrr":
    productTypeBase = 800
elif modelName == "ecmwf-hres":
    productTypeBase = 1000
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
    validTime = initDateTime + timedelta(hours=fhour)
    if "ecmwf" in modelName:
        noticeStr = "Copyright © "+initDateTime.strftime("%Y")+" European Centre for Medium-Range Weather Forecasts (ECMWF)\nhttps://www.ecmwf.int/"
    else:
        noticeStr = None
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n2m Temp, 10m Winds, MSLP", validTime, fhour=fhour, notice=noticeStr, plotHandle=contourmap, cbticks=np.sort(np.append(np.arange(-40, 120, 10), 32)), tickhighlight=32, cbextend="both", colorbarLabel="Temperature (°F)")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output/products/"+modelName+"/sfcTWndMSLP/"+runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig.savefig(staticSavePath+"/f"+str(fhour)+".png", bbox_inches="tight")
    plt.close(fig)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 3
    if hasHelpers:
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)

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
    temperatureColorMap = pltcolors.LinearSegmentedColormap.from_list("temperatureColorMap", all_colors)
    temperatureNorm = pltcolors.TwoSlopeNorm(vcenter=32, vmin=-40, vmax=130)
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
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productTypeBase, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
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
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
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
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
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
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    return windbarbs

def simReflectivityPlot(fileToRead, standaloneFig, ax=None):
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    pathToRead = path.join(inputPath, fileToRead)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    if "refc" in list(modelDataArray.variables):
        prodString = "simrefc"
        simDBZ = modelDataArray.refc
    elif "refd" in list(modelDataArray.variables):
        prodString = "simrefd"
        simDBZ = modelDataArray.refd
    else:
        return
    specR = plt.cm.Spectral_r(np.linspace(0, 1, 200))
    pink = plt.cm.PiYG(np.linspace(0, 0.25, 40))
    purple = plt.cm.PRGn_r(np.linspace(0.75, 1, 40))
    cArr = np.vstack((specR, pink, purple))
    cmap = pltcolors.LinearSegmentedColormap.from_list("cvd-reflectivity", cArr)
    vmin=10
    vmax=80
    dataMask = np.where(np.logical_and(simDBZ.data>=10, simDBZ.data<=80), 0, 1)
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    rdr = ax.pcolormesh(simDBZ.longitude, simDBZ.latitude, np.ma.masked_array(simDBZ, mask=dataMask), cmap=cmap, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree(), zorder=1)

    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, prodString, runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if prodString == "simrefc":
            addon = 8
        elif prodString == "simrefd":
            addon = 80
        productId = productTypeBase + addon
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    return rdr

def updraftHelicityPlot(standaloneFig, ax=None):
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    pathToRead = path.join(inputPath, "udh.grib2")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    udhel = modelDataArray.unknown
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
            lonsToPlot = np.tile(np.array([udhel.longitude.data]), (udhel.data.shape[0], 1))
            latsToPlot = np.tile(udhel.latitude.data, (udhel.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = udhel.longitude
        latsToPlot = udhel.latitude
    ax.contourf(lonsToPlot, latsToPlot, udhel, levels=[50, 999999], cmap="Greys", vmin=0, vmax=100, transform=ccrs.PlateCarree(), zorder=2, transform_first=True, alpha=0.5)
    ax.contour(lonsToPlot, latsToPlot, udhel, levels=[50], colors="black", transform=ccrs.PlateCarree(), zorder=2, transform_first=True, linewidths=0.5)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, "udhelicity", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 9
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)

def staticSimDBZPlot(compOrAGL):
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    sfcWindPlot(False, ax=ax)
    updraftHelicityPlot(False, ax=ax)
    if compOrAGL == "refccomposite":
        addon = 10
        titleStr = "Composite Reflectivity, Updraft Helicity > 50 $m^2 s^{-2}$, 10m Winds"
        rdr = simReflectivityPlot("refc.grib2", False, ax=ax)
    elif compOrAGL == "refdcomposite":
        addon = 81
        titleStr = "1km AGL Reflectivity, Updraft Helicity > 50 $m^2$ $s^{-2}$, 10m Winds"
        rdr = simReflectivityPlot("refd.grib2", False, ax=ax)

    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n"+titleStr, validTime, fhour=fhour, notice=None, plotHandle=rdr, colorbarLabel="Simulated Reflectivity (dBZ)")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output", "products", modelName, compOrAGL, runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
    plt.close(fig)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + addon
    if hasHelpers:
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)


if __name__ == "__main__":
    writeToStatus(str("Plotting init "+str(initDateTime.hour)+"Z f"+str(fhour)+" "+modelName+" "+fieldToPlot))
    inputPath = path.join(basePath, "modelData/"+modelName+"/"+dt.strftime(initDateTime, "%Y%m%d")+"/"+dt.strftime(initDateTime, "%H")+"/"+str(fhour))
    sfcTempPath = path.join(inputPath, "t2m.grib2")
    sfcWindsPath = path.join(inputPath, "sfcwind.grib2")
    sfcPressPath = path.join(inputPath, "sp.grib2")
    windsAtHeightPath = path.join(inputPath, "winds.grib2")
    compositeReflectivityPath = path.join(inputPath, "refc.grib2")
    updraftHelicityPath = path.join(inputPath, "udh.grib2")
    aglReflectivityPath = path.join(inputPath, "refd.grib2")
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
        if fieldToPlot == "refc" and path.exists(compositeReflectivityPath):
            simReflectivityPlot("refc.grib2", True)
        if fieldToPlot == "udh" and path.exists(updraftHelicityPath):
            updraftHelicityPlot(True)
        if fieldToPlot == "refccomposite" and path.exists(compositeReflectivityPath) and path.exists(updraftHelicityPath) and path.exists(sfcWindsPath):
            staticSimDBZPlot("refccomposite")
        if fieldToPlot == "refd" and path.exists(aglReflectivityPath):
            simReflectivityPlot("refd.grib2", True)
        if fieldToPlot == "refdcomposite" and path.exists(aglReflectivityPath) and path.exists(updraftHelicityPath) and path.exists(sfcWindsPath):
            staticSimDBZPlot("refdcomposite")       
