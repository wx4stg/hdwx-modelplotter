#!/usr/bin/env python3
# Forecast model plotting for next-gen HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

from os import fpathconf, path, listdir
from pathlib import Path
import xarray as xr
import multiprocessing as mp
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

def staticSFCTempWindMSLPPlot(pathToRead):
    basePath = path.dirname(path.abspath(__file__))
    modelDataArray = xr.open_dataset(pathToRead)
    modelDataArray = modelDataArray.metpy.parse_cf()
    initTime = dt.strptime(modelDataArray.attrs["_CoordinateModelRunDate"], "%Y-%m-%dT%H:%M:%SZ")
    runPathExt = initTime.strftime("%Y/%m/%d/%H%M")
    staticSavePath = path.join(path.join(basePath, "output/products/hrrr/sfcTwindMSLP/"), runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    if "TMPK_HGHT" in dict(modelDataArray.data_vars).keys():
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
        contourmap = ax.contourf(tempData.longitude, tempData.latitude, tempData, levels=np.arange(-20, 120, 5), cmap="nipy_spectral", vmin=-20, vmax=120, transform=ccrs.PlateCarree(), transform_first=True)
    if "UREL_HGHT" in dict(modelDataArray.data_vars).keys() and "VREL_HGHT" in dict(modelDataArray.data_vars).keys():
        uwind = modelDataArray["UREL_HGHT"]
        uwind = uwind.isel(time=0)
        uwind = uwind.isel(HGHT=0)
        uwind = uwind.metpy.assign_latitude_longitude()
        uwind = uwind.metpy.quantify()
        uwind = uwind.metpy.convert_units("kt")
        vwind = modelDataArray["VREL_HGHT"]
        vwind = vwind.isel(time=0)
        vwind = vwind.isel(HGHT=0)
        vwind = vwind.metpy.assign_latitude_longitude()
        vwind = vwind.metpy.quantify()
        vwind = vwind.metpy.convert_units("kt")
        limit = (slice(None, None, 50), slice(None, None, 50))
        windbarbs = ax.barbs(uwind.longitude.data[limit], uwind.latitude.data[limit], uwind.data[limit], vwind.data[limit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5)
    if "PMSL_NONE" in dict(modelDataArray.data_vars).keys():
        mslpData = modelDataArray["PMSL_NONE"]
        mslpData = mslpData.isel(time=0)
        mslpData = mslpData.metpy.assign_latitude_longitude()
        mslpData = mslpData.metpy.quantify()
        mslpData = mslpData.metpy.convert_units("hPa")
        pressContour = ax.contour(mslpData.longitude, mslpData.latitude, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree(), transform_first=True)
        ax.clabel(pressContour, levels=np.arange(800, 1200, 2), inline=True)
    ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    set_size(1920*px, 1080*px, ax=ax)
    ax.set_extent([modelDataArray.attrs["geospatial_lon_min"], modelDataArray.attrs["geospatial_lon_max"], modelDataArray.attrs["geospatial_lat_min"], modelDataArray.attrs["geospatial_lat_max"]])
    cbax = fig.add_axes([ax.get_position().x0,0.075,(ax.get_position().width/3),.02])
    cb = fig.colorbar(contourmap, cax=cbax, orientation="horizontal")
    cbax.set_xlabel("Temperature (°F)")
    tax = fig.add_axes([ax.get_position().x0+cbax.get_position().width+.01,0.045,(ax.get_position().width/3),.05])
    validTime = str(tempData["time"].data)
    validTime = validTime.split(".")[0]
    validTime = dt.strptime(validTime, "%Y-%m-%dT%H:%M:%S")
    forecastHour = validTime - initTime
    forecastHour = int(forecastHour.seconds / 3600)
    tax.text(0.5, 0.5, initTime.strftime("%H")+"Z HRRR\n2m Temp, 10m Winds, MSLP\nf"+str(forecastHour)+"Valid "+validTime.strftime("%-d %b %Y %H%MZ"), horizontalalignment="center", verticalalignment="center", fontsize=16)
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
    writeJson(303, ["0,0", "0,0"], validTime, forecastHour)
    plt.close(fig)

def tempPlot(pathToRead):
    basePath = path.dirname(path.abspath(__file__))
    modelDataArray = xr.open_dataset(pathToRead)
    if "TMPK_HGHT" not in dict(modelDataArray.data_vars).keys():
        return
    modelDataArray = modelDataArray.metpy.parse_cf()
    initTime = dt.strptime(modelDataArray.attrs["_CoordinateModelRunDate"], "%Y-%m-%dT%H:%M:%SZ")
    runPathExt = initTime.strftime("%Y/%m/%d/%H%M")
    gisTempSavePath = path.join(path.join(basePath, "output/gisproducts/hrrr/sfcT/"), runPathExt)
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
    ax.set_extent([modelDataArray.attrs["geospatial_lon_min"], modelDataArray.attrs["geospatial_lon_max"], modelDataArray.attrs["geospatial_lat_min"], modelDataArray.attrs["geospatial_lat_max"]])
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(path.join(gisTempSavePath, "f"+str(forecastHour)+".png"), transparent=True, bbox_inches=extent)
    gisInfo = [str(modelDataArray.attrs["geospatial_lat_min"])+","+str(modelDataArray.attrs["geospatial_lon_min"]), str(modelDataArray.attrs["geospatial_lat_max"])+","+str(modelDataArray.attrs["geospatial_lon_max"])]
    writeJson(300, gisInfo, validTime, forecastHour)
    plt.close(fig)

def windPlot(pathToRead):
    basePath = path.dirname(path.abspath(__file__))
    modelDataArray = xr.open_dataset(pathToRead)
    if "UREL_HGHT" not in dict(modelDataArray.data_vars).keys():
        return
    if "VREL_HGHT" not in dict(modelDataArray.data_vars).keys():
        return
    modelDataArray = modelDataArray.metpy.parse_cf()
    initTime = dt.strptime(modelDataArray.attrs["_CoordinateModelRunDate"], "%Y-%m-%dT%H:%M:%SZ")
    runPathExt = initTime.strftime("%Y/%m/%d/%H%M")
    gisTempSavePath = path.join(path.join(basePath, "output/gisproducts/hrrr/sfcWnd/"), runPathExt)
    Path(gisTempSavePath).mkdir(parents=True, exist_ok=True)
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    uwind = modelDataArray["UREL_HGHT"]
    uwind = uwind.isel(time=0)
    uwind = uwind.isel(HGHT=0)
    uwind = uwind.metpy.assign_latitude_longitude()
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("kt")
    vwind = modelDataArray["VREL_HGHT"]
    vwind = vwind.isel(time=0)
    vwind = vwind.isel(HGHT=0)
    vwind = vwind.metpy.assign_latitude_longitude()
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("kt")
    limit = (slice(None, None, 50), slice(None, None, 50))
    windbarbs = ax.barbs(uwind.longitude.data[limit], uwind.latitude.data[limit], uwind.data[limit], vwind.data[limit], pivot='middle', color='black', transform=ccrs.PlateCarree(), length=5)

    ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    set_size(1920*px, 1080*px, ax=ax)
    ax.set_extent([modelDataArray.attrs["geospatial_lon_min"], modelDataArray.attrs["geospatial_lon_max"], modelDataArray.attrs["geospatial_lat_min"], modelDataArray.attrs["geospatial_lat_max"]])
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    validTime = str(vwind["time"].data)
    validTime = validTime.split(".")[0]
    validTime = dt.strptime(validTime, "%Y-%m-%dT%H:%M:%S")
    forecastHour = validTime - initTime
    forecastHour = int(forecastHour.seconds / 3600)
    fig.savefig(path.join(gisTempSavePath, "f"+str(forecastHour)+".png"), transparent=True, bbox_inches=extent)
    gisInfo = [str(modelDataArray.attrs["geospatial_lat_min"])+","+str(modelDataArray.attrs["geospatial_lon_min"]), str(modelDataArray.attrs["geospatial_lat_max"])+","+str(modelDataArray.attrs["geospatial_lon_max"])]
    writeJson(302, gisInfo, validTime, forecastHour)
    plt.close(fig)

def mslpPlot(pathToRead):
    basePath = path.dirname(path.abspath(__file__))
    modelDataArray = xr.open_dataset(pathToRead)
    if "PMSL_NONE" not in dict(modelDataArray.data_vars).keys():
        return
    modelDataArray = modelDataArray.metpy.parse_cf()
    initTime = dt.strptime(modelDataArray.attrs["_CoordinateModelRunDate"], "%Y-%m-%dT%H:%M:%SZ")
    runPathExt = initTime.strftime("%Y/%m/%d/%H%M")
    gisTempSavePath = path.join(path.join(basePath, "output/gisproducts/hrrr/mslp/"), runPathExt)
    Path(gisTempSavePath).mkdir(parents=True, exist_ok=True)
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.epsg(3857))
    mslpData = modelDataArray["PMSL_NONE"]
    mslpData = mslpData.isel(time=0)
    mslpData = mslpData.metpy.assign_latitude_longitude()
    mslpData = mslpData.metpy.quantify()
    mslpData = mslpData.metpy.convert_units("hPa")
    pressContour = ax.contour(mslpData.longitude, mslpData.latitude, mslpData, levels=np.arange(800, 1200, 2), colors="black", transform=ccrs.PlateCarree(), transform_first=True)
    ax.clabel(pressContour, levels=np.arange(800, 1200, 2), inline=True)
    ax.add_feature(metpy.plots.USCOUNTIES.with_scale("5m"), edgecolor="gray")
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    set_size(1920*px, 1080*px, ax=ax)
    ax.set_extent([modelDataArray.attrs["geospatial_lon_min"], modelDataArray.attrs["geospatial_lon_max"], modelDataArray.attrs["geospatial_lat_min"], modelDataArray.attrs["geospatial_lat_max"]])
    extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    validTime = str(mslpData["time"].data)
    validTime = validTime.split(".")[0]
    validTime = dt.strptime(validTime, "%Y-%m-%dT%H:%M:%S")
    forecastHour = validTime - initTime
    forecastHour = int(forecastHour.seconds / 3600)
    fig.savefig(path.join(gisTempSavePath, "f"+str(forecastHour)+".png"), transparent=True, bbox_inches=extent)
    gisInfo = [str(modelDataArray.attrs["geospatial_lat_min"])+","+str(modelDataArray.attrs["geospatial_lon_min"]), str(modelDataArray.attrs["geospatial_lat_max"])+","+str(modelDataArray.attrs["geospatial_lon_max"])]
    writeJson(301, gisInfo, validTime, forecastHour)
    plt.close(fig)

if __name__ == "__main__":
    basePath = path.dirname(path.abspath(__file__))
    inputPath = path.join(basePath, "modelData/")
    for file in reversed(sorted(listdir(inputPath))):
        fPath = path.join(inputPath, file)
        try:
            windPlot(fPath)
            mslpPlot(fPath)
            tempPlot(fPath)
            staticSFCTempWindMSLPPlot(fPath)
        except:
            continue


    # inputFiles = [path.join(inputPath, file) for file in sorted(listdir(inputPath))]
    # with mp.Pool(processes=4) as pool:
    #     pool.map(tempPlot, inputFiles)
    #     pool.map(windPlot, inputFiles)