#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 15:05:31 2018

@author: je

Smoothing and Reconstruction of Time Series Vegetation Index with Savitzky and Golday Filter

"""

import glob, os, sys
import rasterio
import gdal, ogr
import numpy as np
import pandas as pd
import datetime as dt
from netCDF4 import Dataset
from itertools import product
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter

def create_TimeSeries (inPath, inVectorFile, outFile):
    """
    - Function to create Vegetation Index Time Serie in NetCDF4 Format 
    - Input should be Folder containing GeoTiff files
    
    """
    # Create Reconstruction Date List
    lstFiles = sorted(glob.glob(inPath+'/*.tif'))
    dicFile = {}
    for File in lstFiles :
#        print (File)
        dicFile.update({os.path.basename(File).split('_')[1]:os.path.basename(File).split('_')[2]})
#    print (dicFile)
    FileName = os.path.basename(lstFiles[0]).split('_')[0]+'_%s_%s_'+\
               os.path.basename(lstFiles[0]).split('_')[3]+'_'+os.path.basename(lstFiles[0]).split('_')[4]
#    print (FileName)
    lstDate =  [Date for Date in dicFile.keys()]  
#    print (lstDate)
    
    # Read Shapefile and Get Extent
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(inVectorFile,gdal.GA_ReadOnly)
    if ds is None:
        print ('Couldn\'t open ShapeFile. Please Check it.')
        sys.exit (1)
    #    print (ds)
    layer = ds.GetLayer()
    extent =  layer.GetExtent()
#    print (extent)
    
    ulx = int(extent[0])
    uly = int(extent[3])
    lrx = int(extent[1])
    lry = int(extent[2])
       
    # Retrieve images profile
    with rasterio.open(lstFiles[0],'r') as ds:
        crs = ds.crs.wkt
        cellsize = ds.transform[1]
        originX = ds.transform[0]
        originY = ds.transform[3]
        lstY = pd.np.arange(lry+0.5*cellsize, uly+0.5*cellsize, cellsize)
        lstX = pd.np.arange(ulx+0.5*cellsize, lrx+0.5*cellsize, cellsize)
          
    # Create Time Series in netCDF4
    
    # Create netcdf file
    ncds = Dataset(outFile ,'w',format='NETCDF4')
#    print (ncds)
    
    # Create Dimensions 
    print ('Creating Dimensions')
#    level = ncds.createDimension("level",None)
    ncds.createDimension("time",len(lstDate))
    ncds.createDimension("X", len(lstX))
    ncds.createDimension("Y", len(lstY))
#    ncds.createDimension("sensor",len(lstDate))#len(lstFiles))
#    print (ncds.dimensions)
    
    # Create Variables
    print ('Create Variables')
    fill_val = -9999.
    crs_var = ncds.createVariable('crs', 'i4')
    crs_var.grid_mapping_name = 'UTM Projection'
    crs_var.proj = crs
    
    y_var = ncds.createVariable('Y', 'f8', ('Y'), fill_value=fill_val)
    y_var.units = 'meters'
    y_var.standard_name = 'Y'

    x_var = ncds.createVariable('X', 'f8', ('X'),fill_value=fill_val)
    x_var.units = 'meters'
    x_var.standard_name = 'X'

    time_var = ncds.createVariable('time', 'l', ('time'),fill_value=fill_val)
    time_var.standard_name = 'time'
    time_var.calendar = 'gregorian'
    
#    sensor_var = ncds.createVariable('sensor', 'S1', ('sensor'),fill_value=fill_val)
#    sensor_var.standard_name = 'Sensor'
#    code_var = ncds.createVariable('code', 'i4', ('latitude', 'longitude'),fill_value=fill_val)

#    outliers_var = ncds.createVariable('outliers', 'f8',('time', 'Y', 'X'),fill_value=fill_val)
#    outliers_var.long_name = 'outliers'

    vi_var = ncds.createVariable('original_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    vi_var.long_name = 'Vegetation Index values'

    savgol_var = ncds.createVariable('savgol_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    savgol_var.long_name = 'hants values'

    combined_var = ncds.createVariable('combined_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    combined_var.long_name = 'combined values'
    
    print ('Load Variables')
    # Load data
    y_var[:] = lstY
    x_var[:] = lstX
    time_var[:] = lstDate
#    sensor_var[:] = df['Sensor'].values
#        
    xOffset = int((ulx- originX) / cellsize)
    yOffset = int((uly - originY) / -cellsize)
    xsize = int((lrx - ulx)/cellsize)  # Raster xsize col
    ysize = int((uly - lry)/cellsize) # Raster ysize lines
    empty_vec = np.empty((len(lstY),len(lstX)))
    empty_vec[:] = -9999.
    
    i = 0
    for Date in lstDate :
        print (Date)
        if (Date in dicFile) : # and dicFile[Date]!='S2')  : # Exclude Sentinel-2 Imagery         
            File = os.path.join(inPath,FileName%(Date,dicFile[Date]))
#            print (File)
            with rasterio.open(File,'r') as ds:
                band = ds.read(1, window=((yOffset, yOffset+ysize),(xOffset, xOffset+xsize)))
    #                    print (band.shape)
                array = np.empty((len(lstY),len(lstX)))
    #                    print (array.shape)
                for j,k in product(range(len(lstY)),range(len(lstX))):
                    array[j,k]=band[len(lstY)-1-j,k]
    #            print (array)
                vi_var [i,:,:] = array
                i+=1
                print (i)
    ncds.close()
    
    print ('########### NetCDF file created #############')
    
    return outFile

def smooth_savgol(ncFile):
    """
    Smooth Time Series in NetCDF4 Format with Savitzky-Golay 
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')
    original_values = ncds.variables['original_values'][:]

    [ztime ,rows, cols] = original_values.shape
    size_st = cols*rows

    values_savgol = np.empty((ztime, rows, cols))

    values_savgol[:] = pd.np.nan

    # Loop
    counter = 1
    print ('Running Savitzky-Golay Filter...')
    for m,n in product(range(rows),range(cols)):
        print ('\t{0}/{1}'.format(counter, size_st))
        y = pd.np.array(original_values[:, m, n])
        yr = savgol_filter(y,window_length=5, polyorder=3)# deriv=0, delta=1.0, axis=-1, mode='interp', cval=0.0)
        values_savgol[:, m, n] = yr
        counter += 1

    ncds.variables['savgol_values'][:] = values_savgol
    # Close netcdf file
    ncds.close()

def check_fit_savgol (ncFile, point):
    """
    Function to check Hants Smoothing
    using Matplotlib in some pixels
    """
    # Read NetCDF File
    ncds = Dataset(ncFile,'r')
    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    
    # Get x and y
    x = point[0]
    y = point[1]
    
#    print (x,y)
    
    # Get row and col index closest to x, y in the netcdf file
    y_closest = Y.flat[pd.np.abs(Y - y).argmin()]
    x_closest = X.flat[pd.np.abs(X - x).argmin()]

    row = pd.np.where(Y == y_closest)[0][0]
    col = pd.np.where(X == x_closest)[0][0]
    
    # Read values
    values_o = ncds.variables['original_values'][:,row,col]
#    print (values_o)
    values_s = ncds.variables['savgol_values'][:,row,col]
#    print (values_s)
    times = [pd.to_datetime(Date, format='%Y%m%d') for Date in ncds.variables['time'][:]]
#    print(times)
    
    plt.figure() 
    
    top = 1.15*max(pd.np.nanmax(values_o), pd.np.nanmax(values_s))
    bottom = min(pd.np.nanmin(values_o), pd.np.nanmin(values_s)) - 0.05
    ylim = [bottom, top]

    # Plot
    plt.plot(times, values_s, 'r-', label='Savitzky-Golay Filter')
    plt.plot(times, values_o, 'b.', label='Original data')

    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc=4)
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: X {0:.2f}, Y {1:.2f} \nFile : {2}'.format(x,y,os.path.basename(ncFile)))
    plt.axes().set_aspect(0.5*(times[-1] - times[0]).days/(ylim[1] - ylim[0]))
    plt.show()

    # Close netcdf file
    ncds.close()

if __name__=="__main__":
    
    inPath = "/home/je/Bureau/Stage/Output/INDICES/NDVI"
    inVectorFile = "/home/je/Bureau/Stage/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"
    outFile = inPath+'/NDVI_savgol.nc'

#    ncFile = create_TimeSeries (inPath,inVectorFile, outFile)
#    ncFile = "/home/je/Bureau/Stage/Output/INDICES/NDVI/NDVI_savgol.nc"
#
#    smooth_savgol (ncFile)


    Points = [(336142.1, 1602889.5), (335329, 1603060), (334582.4, 1600219.6), (334578.5,1600143.5), (336082.9,1603694.0), \
          (336746.3,1603706.7), (336423.9,1603248.0), (336858.4,1603152.2), (337457.6,1603508.4), (337067.1,1603144.3)]
    
    for point in Points :
        check_fit_savgol(ncFile, point)