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
        
        # TRACER sea-breeze convergence thing
        div = mpcalc.divergence(uwind.metpy.convert_units("m/s"), vwind.metpy.convert_units("m/s"))
        ax.scatter([-95.74010, -95.059671], [29.32937, 29.669285], s=50, c="green", marker="s", transform=ccrs.PlateCarree(), zorder=3)
        contourmap = ax.contourf(div.longitude, div.latitude, div, cmap="hot", levels=np.linspace(-.004, 0, 100), transform=ccrs.PlateCarree(), transform_first=True, zorder=1)
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
        if hasHelpers:
            HDWX_helpers.writeJson(basePath, (productTypeBase + 90), initDateTime, "f"+str(fhour)+".png", validTime, ["0,0", "0,0"], 60)
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
