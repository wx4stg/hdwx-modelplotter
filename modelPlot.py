#!/usr/bin/env python3
# Forecast model plotting for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path
from pathlib import Path
import xarray as xr
from metpy import constants
import metpy
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime as dt
from pandas import Timestamp as pdtimestamp
from matplotlib import image as mpimage
from scipy import ndimage

# modelPlot.py <model> <initialization> <fhour>
modelName = sys.argv[1]
initDateTime = dt.strptime(sys.argv[2], "%Y%m%d%H%M")
fhour = int(sys.argv[3])
basePath = path.dirname(path.abspath(__file__))

def writeJson(productID, gisInfo, validTime, fhour):
    productFrameDict = {
        "fhour" : fhour,
        "filename" : "f"+str(fhour)+".png",
        "gisInfo" : gisInfo,
        "valid" : int(validTime.strftime("%Y%m%d%H%M"))
    }

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
    ax = plt.axes(projection=ccrs.epsg(3857))
    set_size(1920*px, 1080*px, ax=ax)
    ax.set_extent([-130, -60, 20, 50])
    contourmap = tempPlot(False, ax=ax)
    validTime = windPlot(False, ax=ax)
    mslpPlot(False, ax=ax)
    ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal")
    cbax.set_xlabel("Temperature (°F)")
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
    staticSavePath = path.join(basePath, "output/products/"+modelName+"/sfcTwindMSLP/"+runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig.savefig(staticSavePath+"/f"+str(fhour)+".png", bbox_inches="tight")
    plt.close(fig)
    gisInfo = ["20,-130", "50,-60"]
    writeJson(304, gisInfo, validTime, fhour)

def tempPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "t2m.grib2")
    modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    modelDataArray = modelDataArray.metpy.parse_cf()
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
    contourmap = ax.contourf(tempData.longitude, tempData.latitude, tempData, levels=np.arange(-20, 120, 5), cmap="nipy_spectral", vmin=-20, vmax=120, transform=ccrs.PlateCarree())
    if standaloneFig:
        ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
        ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
        ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
        set_size(1920*px, 1080*px, ax=ax)
        ax.set_extent([-130, -60, 20, 50])
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    gisInfo = ["20,-130", "50,-60"]
    writeJson(300, gisInfo, validTime, fhour)
    if not standaloneFig:
        return contourmap

def windPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sfcwind.grib2")
    modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    modelDataArray = modelDataArray.metpy.parse_cf()
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
        ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
        ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
        ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
        set_size(1920*px, 1080*px, ax=ax)
        ax.set_extent([-130, -60, 20, 50])
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
        gisInfo = ["20,-130", "50,-60"]
        writeJson(302, gisInfo, validTime, fhour)
    if not standaloneFig:
        return validTime

def mslpPlot(standaloneFig, ax=None):
    pathToRead = path.join(inputPath, "sp.grib2")
    modelDataArray = xr.open_dataset(pathToRead, engine="cfgrib")
    modelDataArray = modelDataArray.metpy.parse_cf()
    runPathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    gisSavePath = path.join(path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/"), runPathExt)
    Path(gisSavePath).mkdir(parents=True, exist_ok=True)
    barometricPressData = modelDataArray["sp"]
    barometricPressData = barometricPressData.metpy.quantify()
    orogData = modelDataArray["orog"]
    orogData = orogData.metpy.quantify()
    tempData = xr.open_dataset(path.join(inputPath, "t2m.grib2"), engine="cfgrib").metpy.parse_cf()
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
    contourmap = ax.contour(barometricPressData.longitude, barometricPressData.latitude, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree())
    if standaloneFig:
        ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
        ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
        ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
        set_size(1920*px, 1080*px, ax=ax)
        ax.set_extent([-130, -60, 20, 50])
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    contourLabels = ax.clabel(contourmap, levels=np.arange(800, 1200, 2), inline=True, fontsize=15)
    [label.set_rotation(0) for label in contourLabels]
    validTime = pdtimestamp(np.datetime64(modelDataArray.valid_time.data)).to_pydatetime()
    gisInfo = ["20,-130", "50,-60"]
    writeJson(301, gisInfo, validTime, fhour)

if __name__ == "__main__":
    inputPath = path.join(basePath, "modelData/"+modelName+"/"+dt.strftime(initDateTime, "%H")+"/"+str(fhour))
    sfcTempPath = path.join(inputPath, "t2m.grib2")
    sfcWindsPath = path.join(inputPath, "sfcwind.grib2")
    sfcPressPath = path.join(inputPath, "sp.grib2")
    canPlotTempWindMSLP = True
    if path.exists(sfcTempPath):
        tempPlot(True)
    else:
        canPlotTempWindMSLP = False
    if path.exists(sfcWindsPath):
        windPlot(True)
    else:
        canPlotTempWindMSLP = False
    if path.exists(sfcPressPath) and path.exists(sfcTempPath):
        mslpPlot(True)
    else:
        canPlotTempWindMSLP = False
    if canPlotTempWindMSLP:
        staticSFCTempWindMSLPPlot()