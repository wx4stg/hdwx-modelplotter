#!/usr/bin/env python3
# Model postprocessing for python-based HDWX
# Created 9 September 2021 by Sam Gardner <stgardner4@tamu.edu>

import sys
from os import path
from pathlib import Path
import xarray as xr
import metpy
from metpy import calc as mpcalc
from cartopy import crs as ccrs
from cartopy import feature as cfeat
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from datetime import datetime as dt, timedelta
from matplotlib import colors as pltcolors
import warnings
import cfgrib

# modelPlot.py
axExtent = [-130, -60, 20, 50]
basePath = path.dirname(path.abspath(__file__))

hasHelpers = False
if path.exists(path.join(basePath, "HDWX_helpers.py")):
    import HDWX_helpers
    hasHelpers = True


def set_size(w,h, ax=None):
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


def staticSFCTempWindMSLPPlot(datasets, modelName, productTypeBase, initDateTime, fhour, runPathExt):
    # Create figure
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())

    # Add artists
    # modelDataArray=datasets['2m'], standaloneFig=True, **model_kwargs
    contourmap = sfcTempPlot(modelDataArray=datasets["2m"], standaloneFig=False, modelName=modelName, ax=ax)
    sfcWindPlot(modelDataArray=datasets["10m"], standaloneFig=False, modelName=modelName, ax=ax)
    if modelName == "ecmwf-hres":
        mslpPlot(modelDataArray=datasets["sfc"], standaloneFig=False, modelName=modelName, ax=ax)
    else:
        mslpPlot(modelDataArray=datasets, standaloneFig=False, modelName=modelName, ax=ax)

    
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n2m Temp, 10m Winds, MSLP", validTime, fhour=fhour, notice=None,
                                plotHandle=contourmap, cbticks=np.sort(np.append(np.arange(-40, 120, 10), 32)), tickhighlight=32, cbextend="both", colorbarLabel="Temperature (°F)")
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


def sfcTempPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
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
        gisSavePath = path.join(basePath, "output/gisproducts/"+modelName+"/sfcT/", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def sfcWindPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
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
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, "sfcWnd", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def mslpPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if modelName == "ecmwf-hres":
        mslpData = modelDataArray["msl"]
        mslpData = mslpData.metpy.quantify()
        mslpData = mslpData.metpy.convert_units("hPa")
    else:
        barometricPressData = modelDataArray["sfc"]["sp"]
        barometricPressData = barometricPressData.metpy.quantify()
        orogData = modelDataArray["sfc"]["orog"]
        orogData = orogData.metpy.quantify()
        tempData = modelDataArray["2m"]["t2m"]
        tempData = tempData.metpy.quantify()
        # I tried using mpcalc altimeter->mslp function here, but it ended up doing nothing and I don't feel like figuring out why
        # Therefore I implemented the same equation manually...
        from metpy import constants
        mslpData = barometricPressData * np.exp(orogData*constants.earth_gravity/(constants.dry_air_gas_constant*tempData))
        mslpData = mslpData.metpy.quantify()
        mslpData = mslpData.metpy.convert_units("hPa")
    # Unidata says smoothing MSLP "a little" is... well they didn't comment on why, they just did it, and it makes the rocky mtns less noisy...
    # https://unidata.github.io/python-gallery/examples/MSLP_temp_winds.html
    from scipy import ndimage
    smoothVal = 5 if modelName in ["namnest", "hrrr"] else 3
    mslpData.data = ndimage.gaussian_filter(mslpData.data, smoothVal)
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
        gisSavePath = path.join(basePath, "output/gisproducts/"+modelName+"/sfcMSLP/", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def windsAtHeightPlot(modelDataArray, pressureLevel, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if modelDataArray.isobaricInhPa.data.shape != ():
        modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
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
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"wind", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def simReflectivityPlot(modelDataArray, standaloneFig, modelName=None, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if "refc" in list(modelDataArray.variables):
        prodString = "simrefc"
        simDBZ = modelDataArray.refc
    elif "refd" in list(modelDataArray.variables):
        prodString = "simrefd"
        simDBZ = modelDataArray.refd
    else:
        return
    import pyart
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


def updraftHelicityPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
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


def staticSimDBZPlot(datasets, modelName, compOrAGL, productTypeBase, initDateTime, fhour, runPathExt):
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    sfcWindPlot(modelDataArray=datasets["10m"], standaloneFig=False, modelName=modelName, ax=ax)
    plotted_udh = False
    if '5km' in datasets.keys() and fhour != 0:
        updraftHelicityPlot(modelDataArray=datasets['5km'], standaloneFig=False, modelName=modelName, ax=ax)
        plotted_udh = True
    if compOrAGL == "refccomposite":
        addon = 10
        if plotted_udh:
            titleStr = "Composite Reflectivity, Updraft Helicity > 50 $m^2 s^{-2}$, 10m Winds"
        else:
            titleStr = "Composite Reflectivity, 10m Winds"
        rdr = simReflectivityPlot(modelDataArray=datasets['atmosphere'], standaloneFig=False, modelName=modelName, ax=ax)
    elif compOrAGL == "refdcomposite":
        addon = 81
        if plotted_udh:
            titleStr = "1km AGL Reflectivity, Updraft Helicity > 50 $m^2 s^{-2}$, 10m Winds"
        else:
            titleStr = "1km AGL Reflectivity, 10m Winds"
        rdr = simReflectivityPlot(modelDataArray=datasets['1km'], standaloneFig=False, modelName=modelName, ax=ax)

    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n"+titleStr, validTime, fhour=fhour, notice=None, plotHandle=rdr, colorbarLabel="Simulated Reflectivity (dBZ)")
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


def heightsPlot(modelDataArray, pressureLevel, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
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
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"hgt", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def tempsPlot(modelDataArray, pressureLevel, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
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
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"temps", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def rhPlot(modelDataArray, pressureLevel, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if pressureLevel not in modelDataArray.isobaricInhPa.values:
        return
    modelDataArray = modelDataArray.sel(isobaricInhPa=pressureLevel)
    if "ecmwf" in modelName:
        modelDataArray = modelDataArray.sel(latitude=slice(54.5, 14.2), longitude=slice(-144.5, -44.5))
    modelDataArray = modelDataArray.metpy.quantify()
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
    contourmap = ax.contourf(lonsToPlot, latsToPlot, rhData.data.clip(0,100), levels=levelsToContour,  cmap=rhcm, vmin=60, vmax=100, transform=ccrs.PlateCarree(), transform_first=True, zorder=1)
    if standaloneFig:
        gisSavePath = path.join(basePath, "output", "gisproducts", modelName, str(pressureLevel)+"rh", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
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


def vort500Plot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    windsAtHeightPlot(modelDataArray=modelDataArray, pressureLevel=500, standaloneFig=False, modelName=modelName, ax=ax)
    heightsPlot(modelDataArray=modelDataArray, pressureLevel=500, standaloneFig=False, modelName=modelName, ax=ax)
    modelDataArray = modelDataArray.sel(isobaricInhPa=500)
    if modelDataArray.u.attrs["GRIB_gridType"] == "regular_ll":
        modelDataArray = modelDataArray.metpy.assign_crs({"grid_mapping_name":"latitude_longitude"})
    elif modelDataArray.u.attrs["GRIB_gridType"] == "lambert":
        with warnings.catch_warnings():
            warnings.simplefilter(action='ignore', category=UserWarning)
            modelDataArray = modelDataArray.metpy.assign_crs({
                    "semi_major_axis": 6371200.0,
                    "semi_minor_axis": 6371200.0,
                    "grid_mapping_name": "lambert_conformal_conic",
                    "standard_parallel": [
                        modelDataArray.u.attrs["GRIB_Latin1InDegrees"],
                    modelDataArray.u.attrs["GRIB_Latin2InDegrees"]
                    ],
                    "latitude_of_projection_origin": modelDataArray.u.attrs["GRIB_LaDInDegrees"],
                    "longitude_of_central_meridian": modelDataArray.u.attrs["GRIB_LoVInDegrees"],
                }).metpy.assign_y_x()
    if modelName == "namnest":
        modelDataArray = modelDataArray.coarsen(y=5, x=5, boundary="trim").mean()
    uwind = modelDataArray.u
    uwind = uwind.metpy.quantify()
    uwind = uwind.metpy.convert_units("m/s")
    vwind = modelDataArray.v
    vwind = vwind.metpy.quantify()
    vwind = vwind.metpy.convert_units("m/s")
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
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


def jetIsotachsPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if 'gh' in modelDataArray.data_vars:
        heightsPlot(modelDataArray=modelDataArray, pressureLevel=250, standaloneFig=False, modelName=modelName, ax=ax)
    windsAtHeightPlot(modelDataArray=modelDataArray, pressureLevel=250, standaloneFig=False, modelName=modelName, ax=ax)
    if modelDataArray.isobaricInhPa.data.shape != ():
        modelDataArray = modelDataArray.sel(isobaricInhPa=250)
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


def temps850Plot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    heightsPlot(modelDataArray=modelDataArray, pressureLevel=850, standaloneFig=False, modelName=modelName, ax=ax)
    windsAtHeightPlot(modelDataArray=modelDataArray, pressureLevel=850, standaloneFig=False, modelName=modelName, ax=ax)
    tempsHandle = tempsPlot(modelDataArray=modelDataArray, pressureLevel=850, standaloneFig=False, modelName=modelName, ax=ax)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    if standaloneFig:
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n 850 hPa Heights, Winds, Temperatures", validTime, fhour=fhour, notice=None, plotHandle=tempsHandle, colorbarLabel="°C")
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


def rh700Plot(datasets, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.LambertConformal())
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "ecmwf-hres":
        mslpPlot(modelDataArray=datasets["sfc"], standaloneFig=False, modelName=modelName, ax=ax)
    else:
        mslpPlot(modelDataArray=datasets, standaloneFig=False, modelName=modelName, ax=ax)
    rhHandle = rhPlot(modelDataArray=datasets['hPa'], pressureLevel=700, standaloneFig=False, modelName=modelName, ax=ax)

    oneThousandHeights = datasets['hPa'].gh.sel(isobaricInhPa=1000).metpy.quantify()
    fiveHundredHeights = datasets['hPa'].gh.sel(isobaricInhPa=500).metpy.quantify()
    thicknessData = fiveHundredHeights - oneThousandHeights
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([thicknessData.longitude.data]), (thicknessData.shape[0], 1))
        latsToPlot = np.tile(thicknessData.latitude.data, (thicknessData.shape[1], 1)).transpose()
    else:
        lonsToPlot = thicknessData.longitude
        latsToPlot = thicknessData.latitude
    coldLevels = np.arange(5340, np.nanmin(thicknessData.data.magnitude)-.01, -60)[::-1]
    hotLevels = np.arange(5460, np.nanmax(thicknessData.data.magnitude)+.01, 60)
    coldContours = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="blue", levels=coldLevels, linewidths=0.5, linestyles="dashdot", transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
    criticalContour = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="red", levels=[5400], linewidths=2, transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
    hotContours = ax.contour(lonsToPlot, latsToPlot, thicknessData, colors="red", levels=hotLevels, linewidths=0.5, linestyles="dashdot", transform=ccrs.PlateCarree(), transform_first=True, zorder=2)
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


def fourPanelPlot(datasets, modelName, initDateTime, fhour, runPathExt, productTypeBase):
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
    ax1, vortHandle = vort500Plot(modelDataArray=datasets['hPa'], standaloneFig=False, modelName=modelName, ax=ax1)
    ax1.set_title("500 hPa Heights, Winds, Absolute Vorticity")
    if 'upper_air' in datasets.keys():
        ax2, jetHandle = jetIsotachsPlot(modelDataArray=datasets['upper_air'], standaloneFig=False, modelName=modelName, ax=ax2)
    else:
        ax2, jetHandle = jetIsotachsPlot(modelDataArray=datasets['hPa'], standaloneFig=False, modelName=modelName, ax=ax2)
    ax2.set_title("250 hPa Heights, Winds, Isotachs")
    ax3, tempHandle = temps850Plot(modelDataArray=datasets['hPa'], standaloneFig=False, modelName=modelName, ax=ax3)
    ax3.set_title("850 hPa Heights, Winds, Temperatures")
    ax4, rhHandle = rh700Plot(datasets=datasets, standaloneFig=False, modelName=modelName, ax=ax4)
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


def staticDewWindMSLPPlot(datasets, modelName, productTypeBase, initDateTime, fhour, runPathExt):
    fig = plt.figure()
    px = 1/plt.rcParams["figure.dpi"]
    fig.set_size_inches(1920*px, 1080*px)
    ax = plt.axes(projection=ccrs.LambertConformal())
    ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    contourmap = sfcDewPlot(modelDataArray=datasets['2m'], standaloneFig=False, modelName=modelName, ax=ax)
    sfcWindPlot(modelDataArray=datasets["10m"], standaloneFig=False, modelName=modelName, ax=ax)
    if modelName == "ecmwf-hres":
        mslpPlot(modelDataArray=datasets["sfc"], standaloneFig=False, modelName=modelName, ax=ax)
    else:
        mslpPlot(modelDataArray=datasets, standaloneFig=False, modelName=modelName, ax=ax)
    ax.add_feature(cfeat.STATES.with_scale("50m"), linewidth=0.5)
    ax.add_feature(cfeat.COASTLINE.with_scale("50m"), linewidth=0.5)
    validTime = initDateTime + timedelta(hours=fhour)
    if hasHelpers:
        HDWX_helpers.dressImage(fig, ax, initDateTime.strftime("%H")+"Z "+modelName.upper()+"\n2m Dew Point Temp, 10m Winds, MSLP", validTime, fhour=fhour, plotHandle=contourmap, cbextend="both", colorbarLabel="Dew Point (°F)")
    staticSavePath = path.join(basePath, "output/products/"+modelName+"/sfcTdWndMSLP/"+runPathExt)
    Path(staticSavePath).mkdir(parents=True, exist_ok=True)
    gisInfo = ["0,0", "0,0"]
    productId = productTypeBase + 5
    if hasHelpers:
        HDWX_helpers.saveImage(fig, path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
        HDWX_helpers.writeJson(basePath, productId, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
    else:
        fig.savefig(path.join(staticSavePath, "f"+str(fhour)+".png"), bbox_inches="tight")
    plt.close(fig)


def sfcDewPlot(modelDataArray, standaloneFig, modelName, productTypeBase=None, initDateTime=None, fhour=None, runPathExt=None, ax=None):
    dewData = modelDataArray["d2m"]
    dewData = dewData.metpy.quantify()
    dewData = dewData.metpy.convert_units("degF")
    if standaloneFig:
        fig = plt.figure()
        px = 1/plt.rcParams["figure.dpi"]
        fig.set_size_inches(1920*px, 1080*px)
        ax = plt.axes(projection=ccrs.epsg(3857))
        ax.set_extent(axExtent, crs=ccrs.PlateCarree())
    if modelName == "gfs" or modelName == "ecmwf-hres":
        lonsToPlot = np.tile(np.array([dewData.longitude.data]), (dewData.data.shape[0], 1))
        latsToPlot = np.tile(dewData.latitude.data, (dewData.data.shape[1], 1)).transpose()
    else:
        lonsToPlot = dewData.longitude
        latsToPlot = dewData.latitude
    levelsToContour = np.arange(10, 91, 5)
    contourmap = ax.contourf(lonsToPlot, latsToPlot, dewData, levels=levelsToContour,  cmap="BrBG", transform=ccrs.PlateCarree(), transform_first=True, extend="both")
    contourLabels = ax.clabel(contourmap, levels=levelsToContour, inline=True, fontsize=10, colors="black")
    [label.set_rotation(0) for label in contourLabels]
    [label.set_text(label.get_text()+"F") for label in contourLabels]
    if standaloneFig:
        gisSavePath = path.join(basePath, "output/gisproducts/"+modelName+"/sfcTd/", runPathExt)
        Path(gisSavePath).mkdir(parents=True, exist_ok=True)
        set_size(1920*px, 1080*px, ax=ax)
        extent = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
        gisInfo = [str(axExtent[2])+","+str(axExtent[0]), str(axExtent[3])+","+str(axExtent[1])]
        validTime = initDateTime + timedelta(hours=fhour)
        if hasHelpers:
            HDWX_helpers.saveImage(fig, path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
            HDWX_helpers.writeJson(basePath, productTypeBase+4, initDateTime, "f"+str(fhour)+".png", validTime, gisInfo, 60)
        else:
            fig.savefig(path.join(gisSavePath, "f"+str(fhour)+".png"), transparent=True, bbox_inches=extent)
        plt.close(fig)
    return contourmap


def plot_all(datasets_unorganized,  modelName, shouldGIS):
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
        raise ValueError("<model> must be 'gfs', 'nam', 'namnest', 'hrrr', or 'ecmwf-hres'")
    datasets = {}
    for ds in datasets_unorganized:
        if 'ecmwf' in modelName:
            ds = ds.sel(latitude=slice(54.5, 14.2), longitude=slice(-144.5, -44.5))
        if 'heightAboveGround' in ds.coords:
            if ds.heightAboveGround.data == 2:
                datasets['2m'] = ds
            elif ds.heightAboveGround.data == 10:
                datasets['10m'] = ds
            elif ds.heightAboveGround.data == 1000:
                datasets['1km'] = ds
            else:
                raise ValueError(f'Unrecognized heightAboveGround: {ds.heightAboveGround.data}')
        elif 'heightAboveGroundLayer' in ds.coords:
            if ds.heightAboveGroundLayer.data == 5000:
                datasets['5km'] = ds
            else:
                raise ValueError(f'Unrecognized heightAboveGroundLayer: {ds.heightAboveGroundLayer.data}')
        elif 'atmosphere' in ds.coords:
            datasets['atmosphere'] = ds
        elif 'atmosphereSingleLayer' in ds.coords:
            datasets['atmosphere'] = ds
        elif 'isobaricInhPa' in ds.coords:
            if np.all(ds.isobaricInhPa.data == 250):
                datasets['upper_air'] = ds
            else:
                datasets['hPa'] = ds
        elif 'surface' in ds.coords or 'meanSea' in ds.coords:
            datasets['sfc'] = ds
        else:
            raise ValueError(f'Unrecognized coordinate in model grib: {ds.coords}')
    initDateTime = ds.time.data.astype('datetime64[s]').astype(dt).item()
    pathExt = initDateTime.strftime("%Y/%m/%d/%H%M")
    fhour = int(ds.step.data/3.6e12)
    model_kwargs = {'modelName' : modelName, 'productTypeBase' : productTypeBase, 'initDateTime' : initDateTime,
                    'fhour' : fhour, 'runPathExt': pathExt}
    print(str("Plotting init "+str(initDateTime.hour)+"Z f"+str(fhour)+" "+modelName))

    if '2m' in datasets.keys() and shouldGIS:
        sfcTempPlot(modelDataArray=datasets['2m'], standaloneFig=True, **model_kwargs)
        sfcDewPlot(modelDataArray=datasets['2m'], standaloneFig=True, **model_kwargs)
    
    if '10m' in datasets.keys() and shouldGIS:
        sfcWindPlot(modelDataArray=datasets['10m'], standaloneFig=True, **model_kwargs)
    
    if ((modelName == "ecmwf-hres" and 'sfc' in datasets.keys()) or ("sfc" in datasets.keys() and "2m" in datasets.keys())) and shouldGIS:
        if modelName == "ecmwf-hres":
            mslpPlot(modelDataArray=datasets['sfc'], standaloneFig=True, **model_kwargs)
        else:
            mslpPlot(modelDataArray=datasets, standaloneFig=True, **model_kwargs)
    
    if '2m' in datasets.keys() and '10m' in datasets.keys() and 'sfc' in datasets.keys():
        staticSFCTempWindMSLPPlot(datasets=datasets, modelName=modelName, productTypeBase=productTypeBase, initDateTime=initDateTime, fhour=fhour, runPathExt=pathExt)
        staticDewWindMSLPPlot(datasets=datasets, modelName=modelName, productTypeBase=productTypeBase, initDateTime=initDateTime, fhour=fhour, runPathExt=pathExt)

    if 'hPa' in datasets.keys():
        vort500Plot(modelDataArray=datasets['hPa'], standaloneFig=True, **model_kwargs)
        temps850Plot(modelDataArray=datasets['hPa'], standaloneFig=True, **model_kwargs)
        rh700Plot(datasets=datasets, standaloneFig=True, **model_kwargs)
        fourPanelPlot(datasets=datasets, modelName=modelName, initDateTime=initDateTime, fhour=fhour, runPathExt=pathExt, productTypeBase=productTypeBase)
        if 'upper_air' in datasets.keys():
            jetIsotachsPlot(modelDataArray=datasets['upper_air'], standaloneFig=True, **model_kwargs)
        else:
            jetIsotachsPlot(modelDataArray=datasets['hPa'], standaloneFig=True, **model_kwargs)
        if shouldGIS:
            [tempsPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [850]]
            [rhPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [700]]
            if 'upper_air' in datasets.keys():
                windsAtHeightPlot(modelDataArray=datasets['upper_air'], pressureLevel=250, standaloneFig=True, **model_kwargs)
                [windsAtHeightPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [500, 850]]
                [heightsPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [500, 850]]
            else:
                [windsAtHeightPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [250, 500, 850]]
                [heightsPlot(modelDataArray=datasets['hPa'], pressureLevel=pressSfc, standaloneFig=True, **model_kwargs) for pressSfc in [250, 500, 850]]

    if 'atmosphere' in datasets.keys():
        staticSimDBZPlot(datasets=datasets, modelName=modelName, compOrAGL='refccomposite', productTypeBase=productTypeBase, initDateTime=initDateTime, fhour=fhour, runPathExt=pathExt)
        if shouldGIS:
            simReflectivityPlot(modelDataArray=datasets['atmosphere'], standaloneFig=True, **model_kwargs)
    if '5km' in datasets.keys() and fhour != 0 and shouldGIS:
        updraftHelicityPlot(modelDataArray=datasets['5km'], standaloneFig=True, **model_kwargs)

    if '1km' in datasets.keys():
        staticSimDBZPlot(datasets=datasets, modelName=modelName, compOrAGL='refdcomposite', productTypeBase=productTypeBase, initDateTime=initDateTime, fhour=fhour, runPathExt=pathExt)
        if shouldGIS:
            simReflectivityPlot(modelDataArray=datasets['1km'], standaloneFig=True, **model_kwargs)
