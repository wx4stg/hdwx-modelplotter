#!/usr/bin/env python3
# Forecast model plotting for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import path, listdir
from pathlib import Path
import xarray as xr
from metpy.units import units
import metpy
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt, rcParams
import numpy as np
from datetime import datetime as dt
from matplotlib import image as mpimage

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


def tempPlot(modelDataArray):
    initTime = dt.strptime(modelDataArray.attrs["_CoordinateModelRunDate"], "%Y-%m-%dT%H:%M:%SZ")
    runPathExt = initTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(path.join(basePath, "output/products/hrrr/sfcTwindMSLP/"), runPathExt)
    gisTempSavePath = path.join(path.join(basePath, "output/gisproducts/hrrr/sfcT/"), runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    Path(gisTempSavePath).mkdir(parents=True, exist_ok=True)
    tempData = modelDataArray["TMPK_HGHT"]
    tempData = tempData.isel(time=0)
    if "HGHT1" in tempData.dims:
        tempData = tempData.isel(HGHT1=0)
    elif "HGHT2" in tempData.dims:
        tempData = tempData.isel(HGHT2=0)
    else:
        print("Height variable not found, this will probably fail shortly...")
        print(tempData.dims)
    tempData = tempData.metpy.assign_latitude_longitude()
    tempData = tempData.metpy.quantify()
    tempData = tempData.metpy.convert_units("degF")
    validTime = str(tempData["time"].data)
    validTime = validTime.split(".")[0]
    validTime = dt.strptime(validTime, "%Y-%m-%dT%H:%M:%S")
    forecastHour = validTime - initTime
    forecastHour = int(forecastHour.seconds / 3600)
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    contourmap = ax.contourf(tempData.longitude, tempData.latitude, tempData, levels=np.arange(-20, 120, 5), cmap="nipy_spectral", vmin=-20, vmax=120, transform=ccrs.PlateCarree(), transform_first=True)
    ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    set_size(1920*px, 1080*px, ax=ax)
    ax.set_extent([modelData.attrs["geospatial_lon_min"], modelData.attrs["geospatial_lon_max"], modelData.attrs["geospatial_lat_min"], modelData.attrs["geospatial_lat_max"]])
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(path.join(gisTempSavePath, "f"+str(forecastHour)+".png"), transparent=True, bbox_inches=extent)
    gisInfo = [str(modelData.attrs["geospatial_lat_min"])+","+str(modelData.attrs["geospatial_lon_min"]), str(modelData.attrs["geospatial_lat_max"])+","+str(modelData.attrs["geospatial_lon_max"])]
    writeJson(300, gisInfo, validTime, forecastHour)
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal")
    cbax.set_xlabel("Temperature (°F)")
    tax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+.01,0.045,(ax.get_position().width/3),.05])
    timeStr = validTime.strftime("Valid %-d %b %Y %H%MZ")
    tax.text(0.5, 0.5, "HRRR 2m Temperature\n"+timeStr, horizontalalignment="center", verticalalignment="center", fontsize=16)
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
    fig.savefig(path.join(staticSavePath, "f"+str(forecastHour)+".png"), bbox_inches="tight")

if __name__ == "__main__":
    basePath = path.dirname(path.abspath(__file__))
    inputPath = path.join(basePath, "modelData/")
    inputFiles = listdir(inputPath)
    for file in inputFiles:
        modelData = xr.open_dataset(path.join(inputPath, file))
        modelData = modelData.metpy.parse_cf()
        tempPlot(modelData)