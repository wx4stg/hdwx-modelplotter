#!/usr/bin/env python3
# Model postprocessing for python-based HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path, listdir, remove
from pathlib import Path
import xarray as xr
from metpy import constants
import metpy
from metpy import calc as mpcalc
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from datetime import datetime as dt, timedelta
from matplotlib import colors as pltcolors
from scipy import ndimage
import pyart
import warnings

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
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 3
    if hasHelpers:
        HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    else:
        fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
    plt.close(fig)

def sfcTempPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "t2m.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output/gisproducts/"+modelName+"/sfcT/", runPathExt)
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
    contourmap = ax.contourf(lonsToPlot, latsToPlot, tempData, levels=levelsToContour,  cmap=temperatureColorMap, norm=temperatureNorm, transform=ccrs.PlateCarree(), transform_first=True, extend="both")
    ax.contour(lonsToPlot, latsToPlot, tempData, levels=[32], colors="red", transform=ccrs.PlateCarree(), transform_first=True, linewidths=1)
    contourLabels = ax.clabel(contourmap, levels=np.arange(-40, 120, 20), inline=True, fontsize=10, colors="black")
    [label.set_rotation(0) for label in contourLabels]
    [label.set_text(label.get_text()+"F") for label in contourLabels]
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productTypeBase, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
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
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 1
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    return windbarbs

def mslpPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sp.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/", runPathExt)
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
    levelsToContour = np.arange((np.nanmin(mslpData.data) // 2) * 2, np.nanmax(mslpData.data)+2, 2)
    contourmap = ax.contour(lonsToPlot, latsToPlot, mslpData, levels=levelsToContour, colors="black", transform=ccrs.PlateCarree(), transform_first=True, linewidths=0.5)
    contourLabels = ax.clabel(contourmap, levels=levelsToContour, inline=True, fontsize=10)
    [label.set_rotation(0) for label in contourLabels]
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 2
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
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
        spatialLimit = slice(None, None, 10)
        dataLimit = (slice(None, None, 10), slice(None, None, 10))
    elif modelName == "nam":
        spatialLimit = (slice(None, None, 20), slice(None, None, 20))
        dataLimit = (slice(None, None, 20), slice(None, None, 20))
    elif modelName == "namnest" or modelName == "hrrr":
        spatialLimit = (slice(None, None, 80), slice(None, None, 80))
        dataLimit = (slice(None, None, 80), slice(None, None, 80))
    windbarbs = ax.barbs(uwind.longitude.data[spatialLimit], uwind.latitude.data[spatialLimit], uwind.data[dataLimit], vwind.data[dataLimit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5, linewidth=0.5, zorder=2)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
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
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
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
    cmap = "pyart_ChaseSpectral"
    vmin=-10
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
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if prodString == "simrefc":
            addon = 8
        elif prodString == "simrefd":
            addon = 80
        productId = productTypeBase + addon
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
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
    if np.nanmax(udhel.data) > 50:
        ax.contourf(lonsToPlot, latsToPlot, udhel, levels=[50, 999999], cmap="Greys", vmin=0, vmax=100, transform=ccrs.PlateCarree(), zorder=2, transform_first=True, alpha=0.5)
        ax.contour(lonsToPlot, latsToPlot, udhel, levels=[50], colors="black", transform=ccrs.PlateCarree(), zorder=2, transform_first=True, linewidths=0.5)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, "udhelicity", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        productId = productTypeBase + 9
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)

def staticSimDBZPlot(compOrAGL):
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    sfcWindPlot(False, ax=ax)
    if modelName in ["namnest", "hrrr"] and fhour != 0:
        updraftHelicityPlot(False, ax=ax)
    if compOrAGL == "refccomposite":
        addon = 10
        if modelName in ["namnest", "hrrr"]:
            titleStr = "Composite Reflectivity, Updraft Helicity > 50 $m^2 s^{-2}$, 10m Winds"
        else:
            titleStr = "Composite Reflectivity, 10m Winds"
        rdr = simReflectivityPlot("refc.grib2", False, ax=ax)
    elif compOrAGL == "refdcomposite":
        addon = 81
        if modelName in ["namnest", "hrrr"]:
            titleStr = "1km AGL Reflectivity, Updraft Helicity > 50 $m^2 s^{-2}$, 10m Winds"
        else:
            titleStr = "1km AGL Reflectivity, 10m Winds"
        rdr = simReflectivityPlot("refd.grib2", False, ax=ax)

    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n"+titleStr, validTime, fhour=fhour, notice=None, plotHandle=rdr, colorbarLabel="Simulated Reflectivity (dBZ)")
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output", "products", modelName, compOrAGL, runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + addon
    if hasHelpers:
        HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    else:
        fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
    plt.close(fig)

def heightsPlot(pressureLevel, standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "heights.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"hgt", runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    heightData = modelDataArray.gh
    heightData = heightData.metpy.quantify()
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([heightData.longitude.data]), (heightData.data.shape[0], 1))
        latsToPlot = np.tile(heightData.latitude.data, (heightData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = heightData.longitude
        latsToPlot = heightData.latitude
    if pressureLevel < 400:
        contourLevels = np.arange((np.nanmin(heightData.data.magnitude // 100)*100 - 100), np.nanmax(heightData.data.magnitude)+121, 120)
    else:
        contourLevels = np.arange((np.nanmin(heightData.data.magnitude // 100)*100 - 100), np.nanmax(heightData.data.magnitude)+61, 60)
    contourmap = ax.contour(lonsToPlot, latsToPlot, heightData.data, levels=contourLevels, colors="black", transform=ccrs.PlateCarree(), transform_first=True, linewidths=0.5, zorder=2)
    contourLabels = ax.clabel(contourmap, levels=contourLevels, inline=True, fontsize=10)
    [label.set_rotation(0) for label in contourLabels]
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if pressureLevel == 250:
            addon = 22
        elif pressureLevel == 500:
            addon = 18
        elif pressureLevel == 850:
            addon = 26
        productId = productTypeBase + addon
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    return contourmap

def tempsPlot(pressureLevel, standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "temps.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    if "ecmwf" in modelName:
        modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"temps", runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    tempsData = modelDataArray.t
    tempsData = tempsData.metpy.quantify()
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([tempsData.longitude.data]), (tempsData.data.shape[0], 1))
        latsToPlot = np.tile(tempsData.latitude.data, (tempsData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = tempsData.longitude
        latsToPlot = tempsData.latitude
    frozenColorMap = plt.cm.cool_r(np.linspace(0, 1, 25))
    meltedColorMap = plt.cm.turbo(np.linspace(0.35, 1, 25))
    all_colors = np.vstack((frozenColorMap, meltedColorMap))
    temperatureColorMap = pltcolors.LinearSegmentedColormap.from_list("temperatureColorMap", all_colors)
    temperatureNorm = pltcolors.TwoSlopeNorm(vcenter=0, vmin=-40, vmax=40)
    levelsToContour = np.arange(-40, 41, 5)
    tempsData = tempsData.data.to("degC")
    contourmap = ax.contourf(lonsToPlot, latsToPlot, tempsData, levels=levelsToContour,  cmap=temperatureColorMap, norm=temperatureNorm, transform=ccrs.PlateCarree(), transform_first=True, extend="both")
    ax.contour(lonsToPlot, latsToPlot, tempsData, levels=[0], colors="red", transform=ccrs.PlateCarree(), transform_first=True, linewidths=1)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if pressureLevel == 850:
            addon = 27
        productId = productTypeBase + addon
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    return contourmap

def rhPlot(pressureLevel, standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "rh.grib2")
    [remove(path.join(inputPath, psblIdxFile)) if psblIdxFile.endswith("idx") else None for psblIdxFile in listdir(inputPath)]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    if "ecmwf" in modelName:
        modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
    modelDataArray = modelDataArray.metpy.quantify()
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"rh", runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    if "r" in modelDataArray.variables:
        rhData = modelDataArray.r
    elif "dpt" in modelDataArray.variables and "t" in modelDataArray.variables:
        rhData = metpy.calc.relative_humidity_from_dewpoint(modelDataArray.t, modelDataArray.dpt) * 100
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([rhData.longitude.data]), (rhData.data.shape[0], 1))
        latsToPlot = np.tile(rhData.latitude.data, (rhData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = rhData.longitude
        latsToPlot = rhData.latitude
    levelsToContour = np.arange(60, 110, 10)
    rhArr = np.array([[1, 185/255, 0, 1], [0, 1, 0, 1], [0, 200/255, 0, 1], [0, 139/255, 0, 1]])
    rhcm = pltcolors.LinearSegmentedColormap.from_list("hdwx-humidity", rhArr)
    contourmap = ax.contourf(lonsToPlot, latsToPlot, rhData.clip(0, 100), levels=levelsToContour,  cmap=rhcm, vmin=60, vmax=100, transform=ccrs.PlateCarree(), transform_first=True, zorder=1)
    if standaloneFig:
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        if pressureLevel == 700:
            addon = 29
        productId = productTypeBase + addon
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    return contourmap

def vort500Plot(standaloneFig, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    windsAtHeightPlot(500, False, ax=ax)
    heightsPlot(500, False, ax=ax)
    pathToRead = path.join(inputPath, "winds.grib2")
    modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib").sel(isobaricInhPa=500)
    uwind = modelDataArray.u
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("m/s")
    vwind = modelDataArray.v
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("m/s")
    vortData = mpcalc.absolute_vorticity(uwind, vwind)
    posVortMap = plt.cm.plasma_r(np.linspace(0, 1, 180))
    neutralVortMap = np.array([[1, 1, 1, 1]]*90)
    negVortMap = plt.cm.cool_r(np.linspace(0.4, 1, 90))
    vortArr = np.vstack((negVortMap, neutralVortMap, posVortMap))
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([vortData.longitude.data]), (vortData.data.shape[0], 1))
        latsToPlot = np.tile(vortData.latitude.data, (vortData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = vortData.longitude
        latsToPlot = vortData.latitude
    vortcm = pltcolors.LinearSegmentedColormap.from_list("hdwx-vorticity", vortArr)
    vortmap = ax.contourf(lonsToPlot, latsToPlot, vortData.data, levels=np.arange(-.00009, .00028, .00003), cmap=vortcm, transform=ccrs.PlateCarree(), transform_first=True, extend="both", zorder=1)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    if standaloneFig:
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n 500 hPa Heights, Winds, Absolute Vorticity", validTime, fhour=fhour, notice=None, plotHandle=vortmap, colorbarLabel=r"$\frac{1}{s}$")
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        staticSavePath = path.join(basePath, "output", "products", modelName, "500staticvort", runPathExt)
        Path(staticSavePath).mkdir(parents=True, exist_ok=True)
        gisInfo = ["0,0", "0,0"]
        productId = productTypeBase + 20
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        plt.close(fig)
    return ax, vortmap

def jetIsotachsPlot(standaloneFig, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    heightsPlot(250, False, ax=ax)
    windsAtHeightPlot(250, False, ax=ax)
    pathToRead = path.join(inputPath, "winds.grib2")
    modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib").sel(isobaricInhPa=250)
    uwind = modelDataArray.u
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("kt")
    vwind = modelDataArray.v
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("kt")
    speed = mpcalc.wind_speed(uwind, vwind)
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([speed.longitude.data]), (speed.shape[0], 1))
        latsToPlot = np.tile(speed.latitude.data, (speed.shape[1], 1)).transpose()
    else:
        lonsToPlot = speed.longitude
        latsToPlot = speed.latitude
    jetmap = ax.contourf(lonsToPlot, latsToPlot, speed, cmap="plasma_r", levels=np.arange(70, 171, 20), transform=ccrs.PlateCarree(), transform_first=True, extend="max", zorder=1)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    if standaloneFig:
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n 250 hPa Heights, Winds, Isotachs", validTime, fhour=fhour, notice=None, plotHandle=jetmap, colorbarLabel="kt")
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        staticSavePath = path.join(basePath, "output", "products", modelName, "250staticjet", runPathExt)
        Path(staticSavePath).mkdir(parents=True, exist_ok=True)
        gisInfo = ["0,0", "0,0"]
        productId = productTypeBase + 24
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        plt.close(fig)
    return ax, jetmap

def temps850Plot(standaloneFig, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    heightsPlot(850, False, ax=ax)
    windsAtHeightPlot(850, False, ax=ax)
    tempsHandle = tempsPlot(850, False, ax=ax)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    if standaloneFig:
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n 850 hPa Heights, Winds, Temperatures", validTime, fhour=fhour, notice=None, plotHandle=tempsHandle, colorbarLabel="°C")
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        staticSavePath = path.join(basePath, "output", "products", modelName, "850statictemps", runPathExt)
        Path(staticSavePath).mkdir(parents=True, exist_ok=True)
        gisInfo = ["0,0", "0,0"]
        productId = productTypeBase + 28
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        plt.close(fig)
    return ax, tempsHandle

def rh700Plot(standaloneFig, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    mslpPlot(False, ax=ax)
    rhHandle = rhPlot(700, False, ax=ax)
    pathToRead = path.join(inputPath, "heights.grib2")
    heightsData = xr.open_dataset(pathToRead, engine="cfgrib")
    oneThousandHeights = heightsData.sel(isobaricInhPa=1000).metpy.quantify()
    fiveHundredHeights = heightsData.sel(isobaricInhPa=500).metpy.quantify()
    thicknessData = fiveHundredHeights.gh - oneThousandHeights.gh
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([thicknessData.longitude.data]), (thicknessData.shape[0], 1))
        latsToPlot = np.tile(thicknessData.latitude.data, (thicknessData.shape[1], 1)).transpose()
    else:
        lonsToPlot = thicknessData.longitude
        latsToPlot = thicknessData.latitude
    coldLevels = np.arange(5340, np.nanmin(thicknessData.data.magnitude)-.01, -60)[::-1]
    hotLevels = np.arange(5460, np.nanmax(thicknessData.data.magnitude)+.01, 60)
    coldContours = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="blue", levels=coldLevels, linewidths=0.5, transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
    criticalContour = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="red", levels=[5400], linewidths=2, transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
    hotContours = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="red", levels=hotLevels, linewidths=0.5, transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
    coldLabels = ax.clabel(coldContours, levels=coldLevels, inline=True, fontsize=10)
    [label.set_rotation(0) for label in coldLabels]
    hotLabels = ax.clabel(hotContours, levels=hotLevels, inline=True, fontsize=10)
    [label.set_rotation(0) for label in hotLabels]
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    if standaloneFig:
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n 700 hPa RH, MSLP, 1000->500 hPa Thicknesses", validTime, fhour=fhour, notice=None, plotHandle=rhHandle, colorbarLabel="%")
        runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
        staticSavePath = path.join(basePath, "output", "products", modelName, "700staticrh", runPathExt)
        Path(staticSavePath).mkdir(parents=True, exist_ok=True)
        gisInfo = ["0,0", "0,0"]
        productId = productTypeBase + 31
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
            HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        plt.close(fig)
    return ax, rhHandle

def fourPanelPlot():
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    gs = GridSpec(3, 4, figure=fig, width_ratios=[0.01, 1, 1, 0.01], height_ratios=[1, 1, 0.1])
    gs.update(left=75/1920,right=1-(75/1920), top=.9925, bottom=0.01, wspace=0, hspace=0.02)
    cax1 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], projection=ccrs.LambertConformal())
    ax2 = fig.add_subplot(gs[0, 2], projection=ccrs.LambertConformal())
    cax2 = fig.add_subplot(gs[0, 3])
    cax3 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1], projection=ccrs.LambertConformal())
    ax4 = fig.add_subplot(gs[1, 2], projection=ccrs.LambertConformal())
    cax4 = fig.add_subplot(gs[1, 3])
    lax = fig.add_axes([0,0,1/5,.06])
    lax.set_aspect(2821/11071, anchor="SE")
    tax = fig.add_axes([1/3,0,1/4,.11])
    [ax.set_extent(axExtent, crs=ccrs.PlateCarree()) for ax in [ax1, ax2, ax3, ax4]]
    ax1, vortHandle = vort500Plot(False, ax=ax1)
    ax1.set_title("500 hPa Heights, Winds, Absolute Vorticity")
    ax2, jetHandle = jetIsotachsPlot(False, ax=ax2)
    ax2.set_title("250 hPa Heights, Winds, Isotachs")
    ax3, tempHandle = temps850Plot(False, ax=ax3)
    ax3.set_title("850 hPa Heights, Winds, Temperatures")
    ax4, rhHandle = rh700Plot(False, ax=ax4)
    ax4.set_title("700 hPa RH, MSLP, 1000->500 hPa Thicknesses")

    cb1 = fig.colorbar(vortHandle, cax=cax1, orientation="vertical", label="1/s", format="%.0e")
    newticks = []
    for tick in cb1.get_ticks():
        if np.abs(tick) > 1e-10:
            newticks.append(tick)
        else:
            newticks.append(0)
    cb1.set_ticks(newticks)
    cb1.ax.yaxis.set_ticks_position("left")
    cb1.ax.yaxis.set_label_position("left")
    fig.colorbar(jetHandle, cax=cax2, orientation="vertical", label="kt")
    cb3 = fig.colorbar(tempHandle, cax=cax3, orientation="vertical", label="°C")
    cb3.ax.yaxis.set_ticks_position("left")
    cb3.ax.yaxis.set_label_position("left")
    fig.colorbar(rhHandle, cax=cax4, orientation="vertical", label="%")

    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, None, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\nForecast 4-Panel", validTime, fhour=fhour, notice=None, plotHandle=None, tax=tax, lax=lax)
    tax.set_position([(0.5-tax.get_position().width/2), tax.get_position().y0, tax.get_position().width, tax.get_position().height])
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(basePath, "output", "products", modelName, "4pnl", runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 32
    if hasHelpers:
        HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    else:
        fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
    plt.close(fig)


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
    heightsPath = path.join(inputPath, "heights.grib2")
    tempsPath = path.join(inputPath, "temps.grib2")
    dewpointPath = path.join(inputPath, "dwpt.grib2")
    rhPath = path.join(inputPath, "rh.grib2")
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
        if fieldToPlot == "refccomposite" and path.exists(compositeReflectivityPath) and path.exists(sfcWindsPath):
            if modelName in ["namnest", "hrrr"] and fhour != 0:
                if path.exists(updraftHelicityPath):
                    staticSimDBZPlot("refccomposite")
            else:
                staticSimDBZPlot("refccomposite")
        if fieldToPlot == "refd" and path.exists(aglReflectivityPath):
            simReflectivityPlot("refd.grib2", True)
        if fieldToPlot == "refdcomposite" and path.exists(aglReflectivityPath) and path.exists(sfcWindsPath):
            if modelName in ["namnest", "hrrr"] and fhour != 0:
                if path.exists(updraftHelicityPath) or fhour == 0:
                    staticSimDBZPlot("refdcomposite")
            else:
                staticSimDBZPlot("refdcomposite")
        if fieldToPlot == "heights" and path.exists(heightsPath):
            [heightsPlot(pressSfc, True) for pressSfc in [250, 500, 850]]
        if fieldToPlot == "temps" and path.exists(tempsPath):
            [tempsPlot(pressSfc, True) for pressSfc in [850]]
        if fieldToPlot == "rh" and path.exists(rhPath):
            [rhPlot(pressSfc, True) for pressSfc in [700]]
        if fieldToPlot == "500vort" and path.exists(heightsPath) and path.exists(windsAtHeightPath):
            vort500Plot(True)
        if fieldToPlot == "jetisotachs" and path.exists(heightsPath) and path.exists(windsAtHeightPath):
            jetIsotachsPlot(True)
        if fieldToPlot == "850temps" and path.exists(heightsPath) and path.exists(windsAtHeightPath) and path.exists(tempsPath):
            temps850Plot(True)
        if fieldToPlot == "700rh" and path.exists(sfcPressPath) and path.exists(heightsPath) and path.exists(rhPath):
            rh700Plot(True)
        if fieldToPlot == "4pnl" and path.exists(heightsPath) and path.exists(windsAtHeightPath) and path.exists(tempsPath) and path.exists(sfcPressPath) and path.exists(rhPath):
            fourPanelPlot()
