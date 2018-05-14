#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 09:49:27 2018

@author: je

Smoothing and Reconstruction of Time Series Vegetation Index with HANTS Algorithm

"""

import glob, os, sys
import rasterio
import gdal, ogr
import numpy as np
import pandas as pd
import datetime as dt
import math
from netCDF4 import Dataset
from itertools import product
from copy import deepcopy
import matplotlib.pyplot as plt
#from multiprocessing import Pool
#import time

def create_TimeSeries (inPath, inVectorFile):
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
#    print (lstFiles)
    start = dt.date(int(os.path.basename(lstFiles[0]).split('_')[1][:4]),int(os.path.basename(lstFiles[0]).split('_')[1][4:6]),
                    int(os.path.basename(lstFiles[0]).split('_')[1][6:]))
    
    end = dt.date(int(os.path.basename(lstFiles[-1]).split('_')[1][:4]),int(os.path.basename(lstFiles[-1]).split('_')[1][4:6]),
                    int(os.path.basename(lstFiles[-1]).split('_')[1][6:]))
#    print (start,end)
    
    dtDate = pd.date_range(start,end,freq='D')
    lstDate = [Date.strftime('%Y%m%d') for Date in dtDate]
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
    outFile = inPath+'/NDVI.nc'
    ncds = Dataset(outFile ,'w',format='NETCDF4')
#    print (ncds)
    
    # Create Dimensions 
    print ('Create Dimensions')
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

    outliers_var = ncds.createVariable('outliers', 'f8',('time', 'Y', 'X'),fill_value=fill_val)
    outliers_var.long_name = 'outliers'

    vi_var = ncds.createVariable('original_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    vi_var.long_name = 'Vegetation Index values'

    hants_var = ncds.createVariable('hants_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    hants_var.long_name = 'hants values'

    combined_var = ncds.createVariable('combined_values', 'f8', ('time', 'Y', 'X'), fill_value=fill_val)
    combined_var.long_name = 'combined values'
    
    print ('Load Variables')
    # Load data
    y_var[:] = lstY
    x_var[:] = lstX
    time_var[:] = lstDate
#    sensor_var[:] = df['Sensor'].values
        
    xOffset = int((ulx- originX) / cellsize)
    yOffset = int((uly - originY) / -cellsize)
    xsize = int((lrx - ulx)/cellsize)  # Raster xsize col
    ysize = int((uly - lry)/cellsize) # Raster ysize lines
    empty_vec = np.empty((len(lstY),len(lstX)))
    empty_vec[:] = -9999.
    
    i = 0
    for Date in lstDate :
        print (Date)
#        File = os.path.join(inPath,FileName%Date)
#        if File in lstFiles :
        if Date in dicFile :          
            File = os.path.join(inPath,FileName%(Date,dicFile[Date]))
            print (File)
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
        else :
            vi_var [i,:,:] = empty_vec
            i+=1
            print (i)
                
    ncds.close()
    
    print ('########### NetCDF file created #############')
    
    return outFile
    
def smooth_hants (ncFile, nb, nf, HiLo, low, high, fet, dod, delta, fill_val):
    """
    Smooth Time Series in NetCDF4 Format with Harmonic Analysis of Time series
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')

    time_var = ncds.variables['time'][:]
    original_values = ncds.variables['original_values'][:]

    [ztime ,rows, cols] = original_values.shape
    size_st = cols*rows

    values_hants = np.empty((ztime, rows, cols))
    outliers_hants = np.empty((ztime, rows, cols))

    values_hants[:] = pd.np.nan
    outliers_hants[:] = pd.np.nan

    # Additional parameters
    ni = len(time_var)
    ts = range(ni)

    # Loop
    counter = 1
    print ('Running HANTS...')
    for m,n in product(range(rows),range(cols)):
        print ('\t{0}/{1}'.format(counter, size_st))
    
        y = pd.np.array(original_values[:, m, n])
    
        y[pd.np.isnan(y)] = fill_val
    
        [yr, outliers] = HANTS(ni, nb, nf, y, ts, HiLo,
                               low, high, fet, dod, delta, fill_val)
    
        values_hants[:, m, n] = yr
        outliers_hants[:, m, n] = outliers
    
        counter = counter + 1

    ncds.variables['hants_values'][:] = values_hants
    ncds.variables['outliers'][:] = outliers_hants
    ncds.variables['combined_values'][:] = pd.np.where(outliers_hants,values_hants,original_values)
    # Close netcdf file
    ncds.close()

def check_fit_hants (ncFile, point):
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
    values_h = ncds.variables['hants_values'][:,row,col]
#    print (values_h)
    times = [pd.to_datetime(Date, format='%Y%m%d') for Date in ncds.variables['time'][:]]
#    print(times)
    
    plt.figure() 
    
    top = 1.15*max(pd.np.nanmax(values_o), pd.np.nanmax(values_h))
    bottom = min(pd.np.nanmin(values_o), pd.np.nanmin(values_h)) - 0.05
    ylim = [bottom, top]

    # Plot
    plt.plot(times, values_h, 'r-', label='HANTS')
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
  
def HANTS(ni, nb, nf, y, ts, HiLo, low, high, fet, dod, delta, fill_val):
    """
    This function applies the Harmonic ANalysis of Time Series (HANTS)
    algorithm originally developed by the Netherlands Aerospace Centre (NLR)
    (http://www.nlr.org/space/earth-observation/).

    This python implementation was based on two previous implementations
    available at the following links:
    https://codereview.stackexchange.com/questions/71489/harmonic-analysis-of-time-series-applied-to-arrays
    http://nl.mathworks.com/matlabcentral/fileexchange/38841-matlab-implementation-of-harmonic-analysis-of-time-series--hants-
    """
    # Arrays
    mat = pd.np.zeros((min(2*nf+1, ni), ni))
    # amp = np.zeros((nf + 1, 1))

    # phi = np.zeros((nf+1, 1))
    yr = pd.np.zeros((ni, 1))
    outliers = pd.np.zeros((1, len(y)))

    # Filter
    sHiLo = 0
    if HiLo == 'Hi':
        sHiLo = -1
    elif HiLo == 'Lo':
        sHiLo = 1

    nr = min(2*nf+1, ni)
    noutmax = ni - nr - dod
    # dg = 180.0/math.pi
    mat[0, :] = 1.0

    ang = 2*math.pi*pd.np.arange(nb)/nb
    cs = pd.np.cos(ang)
    sn = pd.np.sin(ang)

    i = pd.np.arange(1, nf+1)
    for j in pd.np.arange(ni):
        index = pd.np.mod(i*ts[j], nb)
        mat[2 * i-1, j] = cs.take(index)
        mat[2 * i, j] = sn.take(index)

    p = pd.np.ones_like(y)
    bool_out = (y < low) | (y > high)
    p[bool_out] = 0
    outliers[bool_out.reshape(1, y.shape[0])] = 1
    nout = pd.np.sum(p == 0)

    if nout > noutmax:
        if pd.np.isclose(y, fill_val).any():
            ready = pd.np.array([True])
            yr = y
            outliers = pd.np.zeros((y.shape[0]), dtype=int)
            outliers[:] = fill_val
        else:
            raise Exception('Not enough data points.')
    else:
        ready = pd.np.zeros((y.shape[0]), dtype=bool)

    nloop = 0
    nloopmax = ni

    while ((not ready.all()) & (nloop < nloopmax)):

        nloop += 1
        za = pd.np.matmul(mat, p*y)

        A = pd.np.matmul(pd.np.matmul(mat, pd.np.diag(p)),
                         pd.np.transpose(mat))
        A = A + pd.np.identity(nr)*delta
        A[0, 0] = A[0, 0] - delta

        zr = pd.np.linalg.solve(A, za)

        yr = pd.np.matmul(pd.np.transpose(mat), zr)
        diffVec = sHiLo*(yr-y)
        err = p*diffVec

        err_ls = list(err)
        err_sort = deepcopy(err)
        err_sort.sort()

        rankVec = [err_ls.index(f) for f in err_sort]

        maxerr = diffVec[rankVec[-1]]
        ready = (maxerr <= fet) | (nout == noutmax)

        if (not ready):
            i = ni - 1
            j = rankVec[i]
            while ((p[j]*diffVec[j] > 0.5*maxerr) & (nout < noutmax)):
                p[j] = 0
                outliers[0, j] = 1
                nout += 1
                i -= 1
                if i == 0:
                    j = 0
                else:
                    j = 1

    return [yr, outliers]

if __name__=="__main__":
    
    inPath = "/home/je/Bureau/Stage/Output/INDICES/NDVI"
    inVectorFile = "/home/je/Bureau/Stage/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628.shp"
#    inPath = "G:/NDVI"
#    inVectorFile = "H:/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628.shp"
#    ncFile = create_TimeSeries (inPath,inVectorFile)
    
#    nb = 365
#    nf = 3
#    HiLo = 'Lo'
#    low = -0.3
#    high = 1
#    fet = 0.05
#    dod = 1
#    delta = 0.25
#    fill_val = -9999.
#    smooth_hants (ncFile, nb, nf, HiLo, low, high, fet, dod, delta, fill_val)
    
    
#    ncFile = "/home/je/Bureau/Stage/Output/INDICES/NDVI/NDVI.nc" 
    ncFile = "/home/je/Bureau/Stage/Output/INDICES/NDVI/NDVI_Planet.nc"

    Points = [(336142.1, 1602889.5), (335329, 1603060), (334582.4, 1600219.6), (334578.5,1600143.5), (336082.9,1603694.0), \
          (336746.3,1603706.7), (336423.9,1603248.0), (336858.4,1603152.2), (337457.6,1603508.4), (337067.1,1603144.3)]
    for point in Points :
        check_fit_hants(ncFile, point)

        