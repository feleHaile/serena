#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 09:46:49 2018

@author: je

Create Time Series from VI Index in NetCDF4 Format
And Smooth values by Hants, Savitzky-Golay and Whittaker-Henderson
"""
import glob, os, sys
import datetime as dt
import rasterio
import gdal, ogr
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from itertools import product
from copy import deepcopy
import math
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.misc import comb
from scipy import sparse
from scipy.sparse.linalg import spsolve

def create_time_series (inPath, inVectorFile):
    """
    - Function to create Vegetation Index Time Series in NetCDF4 Format 
    - 3 Time Series variables are created : PlanetScope - PlanetScope & RapidEye - PlanetScope, RapidEye & Sentinel-2
    - Input should be Folder containing GeoTiff files
    - VI Index Values have 10 000 scale factor
    """
    # Create Date List
    lstFiles = sorted(glob.glob(inPath+'/*.tif'))
    dicFile = {}
    for File in lstFiles :
#        print (File)
        dicFile.update({os.path.basename(File).split('_')[1]:os.path.basename(File).split('_')[2]})
#    print (dicFile)
    FileName = os.path.basename(lstFiles[0]).split('_')[0]+'_%s_%s_'+\
               os.path.basename(lstFiles[0]).split('_')[3]+'_'+os.path.basename(lstFiles[0]).split('_')[4]
#    print (FileName)
    viIndexName = FileName.split('_')[0]
    start = dt.date(int(os.path.basename(lstFiles[0]).split('_')[1][:4]),int(os.path.basename(lstFiles[0]).split('_')[1][4:6]),
                    int(os.path.basename(lstFiles[0]).split('_')[1][6:]))
    
    end = dt.date(int(os.path.basename(lstFiles[-1]).split('_')[1][:4]),int(os.path.basename(lstFiles[-1]).split('_')[1][4:6]),
                    int(os.path.basename(lstFiles[-1]).split('_')[1][6:]))
#    print (start,end)
    
    dtDate = pd.date_range(start,end,freq='D')
    lstDate = [Date.strftime('%Y%m%d') for Date in dtDate]
    
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
    outFile = os.path.join(inPath,viIndexName+'_TIME_SERIES.nc')
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

    vi_var_P = ncds.createVariable('original_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    vi_var_P.long_name = 'PlanetScope %s Index values'%viIndexName
    
    vi_var_PR = ncds.createVariable('original_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    vi_var_PR.long_name = 'PlanetScope & RapidEye %s Index values'%viIndexName
    
    vi_var_PRS = ncds.createVariable('original_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    vi_var_PRS.long_name = 'PlanetScope, RapidEye & Sentinel-2 %s Index values'%viIndexName
    
#    interpolated_var_P = ncds.createVariable('inteprolated_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
#    interpolated_var_P.long_name = 'Interpolated PlanetScope %s Index values'%viIndexName
#    
#    interpolated_var_PR = ncds.createVariable('inteprolated_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
#    interpolated_var_PR.long_name = 'Interpolated PlanetScope & RapidEye %s Index values'%viIndexName
#    
#    interpolated_var_PRS = ncds.createVariable('inteprolated_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
#    interpolated_var_PRS.long_name = 'Interpolated PlanetScope, RapidEye & Sentinel-2 %s Index values'%viIndexName
    
    print ('Load Variables')
    # Load data
    y_var[:] = lstY
    x_var[:] = lstX
    time_var[:] = lstDate
        
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
            print (File)
            with rasterio.open(File,'r') as ds:
                band = ds.read(1, window=((yOffset, yOffset+ysize),(xOffset, xOffset+xsize)))
#                    print (band.shape)
                array = np.empty((len(lstY),len(lstX)))
#                    print (array.shape)
                for j,k in product(range(len(lstY)),range(len(lstX))):
                    array[j,k]=band[len(lstY)-1-j,k]
#               array[pd.np.isnan(array)] = fill_val
                array = np.where(array==np.nan,-9999.,array*10000)
    #            print (array)
                vi_var_PRS [i,:,:] = array
                print ("PRS")
                if (dicFile[Date]!='S2') :
                    vi_var_PR[i,:,:] = array
                    print ("PR")
                if (dicFile[Date]!='S2' and dicFile[Date]!='RapidEye'):
                    vi_var_P[i,:,:] = array
                    print ("P")
                i+=1
#                print (i)
        else :
            vi_var_PRS [i,:,:] =  empty_vec
            vi_var_PR[i,:,:] =  empty_vec
            vi_var_P[i,:,:] =  empty_vec
            i+=1
#            print (i)
    
#    print ("Load Interpolated")
#    counter = 1
#    for m,n in product(range(len(lstY)),range(len(lstX))):
#        print ('%s/%s'%(counter, len(lstY)*len(lstX)))
#        s_P = vi_var_P[:,m,n].data
#        s_P = np.where(s_P==-9999.,np.nan,s_P/10000)
#        
#        s_PR = vi_var_PR[:,m,n].data
#        s_PR = np.where(s_PR==-9999.,np.nan,s_PR/10000)
#        
#        s_PRS = vi_var_PRS[:,m,n].data
#        s_PRS = np.where(s_PRS==-9999.,np.nan,s_PRS/10000)
#        
#        series_P = pd.Series(s_P).interpolate(method='linear')
#        array_P = np.where(series_P.values==np.nan,-9999.,series_P.values*10000)
#        
#        series_PR = pd.Series(s_PR).interpolate(method='linear')
#        array_PR = np.where(series_PR.values==np.nan,-9999.,series_PR.values*10000)
#        
#        series_PRS = pd.Series(s_PRS).interpolate(method='linear')
#        array_PRS = np.where(series_PRS.values==np.nan,-9999.,series_PRS.values*10000)
#        
#        interpolated_var_P [:,m,n] = array_P
#        interpolated_var_PR [:,m,n] = array_PR
#        interpolated_var_PRS [:,m,n] = array_PRS
#        counter+=1
#            
    ncds.close()
    
    print ('NetCDF file created')
    
    return outFile

# =============================================================================
# Smoothing with HANTS
# =============================================================================
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
    
def smooth_hants (ncFile, nb=365, nf=3, HiLo='Lo', low=-0.3, high=1, fet=0.05, dod=1, delta=0.25, fill_val=-9999.):
    """
    Smooth Time Series in NetCDF4 Format with Harmonic Analysis of Time series
    HANTS parameters are :
        nb = 365
        nf = 3
        HiLo = 'Lo'
        low = -0.3
        high = 1
        fet = 0.05
        dod = 1
        delta = 0.25
        fill_val = -9999.
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')
    
    viIndexName = os.path.basename(ncFile).split('_')[0]
    
    # Create HANTS variables
    hants_var_P = ncds.createVariable('hants_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    hants_var_P.long_name = 'PlanetScope %s HANTS values'
    
    hants_var_PR = ncds.createVariable('hants_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    hants_var_PR.long_name = 'PlanetScope & RapidEye %s HANTS values'
    
    hants_var_PRS = ncds.createVariable('hants_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    hants_var_PRS.long_name = 'PlanetScope, RapidEye and Sentinel-2 %s HANTS values'

    combined_var_P = ncds.createVariable('combined_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    combined_var_P.long_name = 'PlanetScope %s combined values'%viIndexName

    combined_var_PR = ncds.createVariable('combined_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    combined_var_PR.long_name = 'PlanetScope & RapidEye %s combined values'%viIndexName
    
    combined_var_PRS = ncds.createVariable('combined_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    combined_var_PRS.long_name = 'PlanetScope, RapidEye and Sentinel-2 %s combined values'%viIndexName
    
    outliers_var_P = ncds.createVariable('outliers_P_%s_values'%viIndexName, 'i2',('time', 'Y', 'X'),fill_value=fill_val)
    outliers_var_P.long_name = 'PlanetScope %s Outliers values '%viIndexName
    
    outliers_var_PR = ncds.createVariable('outliers_PR_%s_values'%viIndexName, 'i2',('time', 'Y', 'X'),fill_value=fill_val)
    outliers_var_PR.long_name = 'PlanetScope & RapidEye %s Outliers values '%viIndexName
    
    outliers_var_PRS = ncds.createVariable('outliers_PRS_%s_values'%viIndexName, 'i2',('time', 'Y', 'X'),fill_value=fill_val)
    outliers_var_PRS.long_name = 'PlanetScope, RapidEye and Sentinel-2 %s Outliers values '%viIndexName
    
    # Get Existing Values
    time_var = ncds.variables['time'][:]
    original_P_values = ncds.variables['original_P_%s_values'%viIndexName][:]/10000
    original_PR_values = ncds.variables['original_PR_%s_values'%viIndexName][:]/10000
    original_PRS_values = ncds.variables['original_PRS_%s_values'%viIndexName][:]/10000

    [ztime ,rows, cols] = original_P_values.shape
    size_st = cols*rows

    values_hants_P = np.empty((ztime, rows, cols))
    values_hants_PR = np.empty((ztime, rows, cols))
    values_hants_PRS = np.empty((ztime, rows, cols))
    outliers_hants_P = np.empty((ztime, rows, cols))
    outliers_hants_PR = np.empty((ztime, rows, cols))
    outliers_hants_PRS = np.empty((ztime, rows, cols))

    values_hants_P[:] = pd.np.nan
    values_hants_PR[:] = pd.np.nan
    values_hants_PRS[:] = pd.np.nan
    outliers_hants_P[:] = pd.np.nan
    outliers_hants_PR[:] = pd.np.nan
    outliers_hants_PRS[:] = pd.np.nan

    # Additional parameters
    ni = len(time_var)
    ts = range(ni)

    # Loop
    counter = 1
    print ('Running HANTS...')
    for m,n in product(range(rows),range(cols)):
        print ('\t{0}/{1}'.format(counter, size_st))
    
        yp = pd.np.array(original_P_values[:, m, n])
        yp[pd.np.isnan(yp)] = fill_val
        [yr_p, outliers_p] = HANTS(ni, nb, nf, yp, ts, HiLo,low, high, fet, dod, delta, fill_val)
        values_hants_P[:, m, n] = yr_p
        outliers_hants_P[:, m, n] = outliers_p
        
        ypr = pd.np.array(original_PR_values[:, m, n])
        ypr[pd.np.isnan(ypr)] = fill_val
        [yr_pr, outliers_pr] = HANTS(ni, nb, nf, ypr, ts, HiLo,low, high, fet, dod, delta, fill_val)
        values_hants_PR[:, m, n] = yr_pr
        outliers_hants_PR[:, m, n] = outliers_pr
        
        yprs = pd.np.array(original_PRS_values[:, m, n])
        yprs[pd.np.isnan(yprs)] = fill_val
        [yr_prs, outliers_prs] = HANTS(ni, nb, nf, yprs, ts, HiLo,low, high, fet, dod, delta, fill_val)
        values_hants_PRS[:, m, n] = yr_prs
        outliers_hants_PRS[:, m, n] = outliers_prs
    
        counter = counter + 1
    
    values_hants_P = values_hants_P * 10000
    values_hants_PR = values_hants_PR * 10000
    values_hants_PRS = values_hants_PRS * 10000
    
    outliers_hants_P = outliers_hants_P * 10000
    outliers_hants_PR = outliers_hants_PR * 10000
    outliers_hants_PRS = outliers_hants_PRS * 10000
    
    ncds.variables['hants_P_%s_values'%viIndexName][:] = values_hants_P
    ncds.variables['hants_PR_%s_values'%viIndexName][:] = values_hants_PR
    ncds.variables['hants_PRS_%s_values'%viIndexName][:] = values_hants_PRS
    ncds.variables['outliers_P_%s_values'%viIndexName][:] = outliers_hants_P
    ncds.variables['outliers_PR_%s_values'%viIndexName][:] = outliers_hants_PR
    ncds.variables['outliers_PRS_%s_values'%viIndexName][:] = outliers_hants_PRS
    ncds.variables['combined_P_%s_values'%viIndexName][:] = pd.np.where(outliers_hants_P,values_hants_P,original_P_values*10000)
    ncds.variables['combined_PR_%s_values'%viIndexName][:] = pd.np.where(outliers_hants_PR,values_hants_PR,original_PR_values*10000)
    ncds.variables['combined_PRS_%s_values'%viIndexName][:] = pd.np.where(outliers_hants_PRS,values_hants_PRS,original_PRS_values*10000)
    # Close netcdf file
    ncds.close()
    
# =============================================================================
# Smoothing with Savitzky-Golay Filter
# =============================================================================

def smooth_savgol(ncFile, window_length=31, polyorder=4):
    """
    Smooth Time Series in NetCDF4 Format with Savitzky-Golay 
    Savitzky-Golay parameters are :
       windows_length :
       polyorder : 
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')
    viIndexName = os.path.basename(ncFile).split('_')[0]
    fill_val = -9999.
    
    # Create Savitzky-Golay variables
    savgol_var_P = ncds.createVariable('savgol_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    savgol_var_P.long_name = 'PlanetScope %s Savitzky-Golay Filter values'

    savgol_var_PR = ncds.createVariable('savgol_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    savgol_var_PR.long_name = 'PlanetScope & RapidEye %s Savitzky-Golay Filter values'

    savgol_var_PRS = ncds.createVariable('savgol_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    savgol_var_PRS.long_name = 'PlanetScope, RapidEye and Sentinel-2 %s Savitzky-Golay Filter values'
    
    # Get Existing Values
    original_P_values = ncds.variables['original_P_%s_values'%viIndexName][:]/10000
    original_PR_values = ncds.variables['original_PR_%s_values'%viIndexName][:]/10000
    original_PRS_values = ncds.variables['original_PRS_%s_values'%viIndexName][:]/10000
    
    # Filtering
    [ztime ,rows, cols] = original_P_values.shape
    size_st = cols*rows

    values_savgol_P = np.empty((ztime, rows, cols))
    values_savgol_P[:] = pd.np.nan
    values_savgol_PR = np.empty((ztime, rows, cols))
    values_savgol_PR[:] = pd.np.nan
    values_savgol_PRS = np.empty((ztime, rows, cols))
    values_savgol_PRS[:] = pd.np.nan
    
    counter = 1
    print ('Running Savitzky-Golay Filter...')
    for m,n in product(range(1),range(1)): # rows cols 
        print ('\t{0}/{1}'.format(counter, size_st))
        
        yp = pd.Series(original_P_values[:, m, n])
        yp_inter = yp.interpolate(method='linear')
        yr_p = savgol_filter(yp_inter,window_length, polyorder)
        values_savgol_P[:, m, n] = yr_p

        ypr = pd.Series(original_PR_values[:, m, n])
        ypr_inter = ypr.interpolate(method='linear')
        yr_pr = savgol_filter(ypr_inter,window_length, polyorder)
        values_savgol_PR[:, m, n] = yr_pr
        
        yprs = pd.Series(original_PRS_values[:, m, n])
        yprs_inter = yprs.interpolate(method='linear')
        yr_prs = savgol_filter(yprs_inter,window_length, polyorder)
        values_savgol_PRS[:, m, n] = yr_prs
        
        counter += 1
        
    ncds.variables['savgol_P_%s_values'%viIndexName][:] = values_savgol_P*10000
    ncds.variables['savgol_PR_%s_values'%viIndexName][:] = values_savgol_PR*10000
    ncds.variables['savgol_PRS_%s_values'%viIndexName][:] = values_savgol_PRS*10000
    
    # Close netcdf file
    ncds.close()

# =============================================================================
# Smoothing with Whittaker-Henderson
# =============================================================================

def whfilter(a, weights, lamb, p=3):
    """
    Generalized Whittaker-Handerson Graduation Method
    Parameters
    ----------
    a : array-like
          The input array, shape (n,)
    weights : array-like or None
          Weights
    lamb : float
          The relative importance between goodness of fit 
          and smoothness (smoothness increases with lamb).
    p : integer, default 3
          The degree of smoothness. We minimize the p-th 
          finite-differences of the graduated data. Examples:
          p=2 Hodrick-Prescott filter;
          p=3 Whittaker-Henderson method;
          Note: moments 0..p-1 will be conserved by graduation
    
    Returns
    -------
    out : array
          The smoothed data
    References
    ----------
    implementation of scikits.statsmodels.tsa.filters.hp_filter.py
    Alicja S. Nocon & William F. Scott (2012): "An extension of the 
       Whittaker-Henderson method of graduation", Scandinavian 
       Actuarial Journal, 2012:1, 70-79
    Whittaker, E. T. (1922). "On a new method of graduation", 
       Proceedings of the Edinburgh Mathematical Society 41,63-75.
    """
    # input data
    a = np.squeeze(a)
    if a.ndim>1: raise ValueError("input array a must be 1d")
    n = a.size

    # weights
    W = np.squeeze(weights) if weights is not None else np.ones(n)
    if np.any(W==0) or not np.all(np.isfinite(W)): 
      raise ValueError("weights must be non-zero and finite.")
    W = sparse.dia_matrix((W, 0), shape=(n,n))

    # set up difference Matrix K, shape (n-p, n)
    # K_ij = k(j-i),  l=j-i
    # k(l) = (-1)^l Binomial(p,l) if 0<=l<=p else 0
    l = np.arange(p+1)
    k = (-1)**l * comb(p,l)       # same as K_0j
    diags  =np.tile(k,(n,1)).T   # side-diagonal K_i,i+l; n-times k(l)
    offsets=np.arange(p+1)        # index of side-diagonals
    K = sparse.dia_matrix((diags,offsets),shape=(n-p,n)) # K_ij

    # solve quadratic optimization problem 
    return spsolve(W+lamb*K.T.dot(K), W.dot(a))

def smooth_whittaker(ncFile, weights=None, lamb=1000):
    """
    Smooth Time Series in NetCDF4 Format with Whittaker-Henderson Smoother
    Whittaker-Henderson parameters are :
       weights :
       lamb : 
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')
    viIndexName = os.path.basename(ncFile).split('_')[0]
    fill_val = -9999.
    
    # Create Whittaker-Henderson variables
    whittaker_var_P = ncds.createVariable('whittaker_P_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    whittaker_var_P.long_name = 'PlanetScope %s Whittaker-Henderson Smoother values'

    whittaker_var_PR = ncds.createVariable('whittaker_PR_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    whittaker_var_PR.long_name = 'PlanetScope & RapidEye %s Whittaker-Henderson Smoother values'

    whittaker_var_PRS = ncds.createVariable('whittaker_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    whittaker_var_PRS.long_name = 'PlanetScope, RapidEye and Sentinel-2 %s Whittaker-Henderson Smoother values'
    
    # Get Existing Values
    original_P_values = ncds.variables['original_P_%s_values'%viIndexName][:]/10000
    original_PR_values = ncds.variables['original_PR_%s_values'%viIndexName][:]/10000
    original_PRS_values = ncds.variables['original_PRS_%s_values'%viIndexName][:]/10000
    
    # Filtering
    [ztime ,rows, cols] = original_P_values.shape
    size_st = cols*rows

    values_whittaker_P = np.empty((ztime, rows, cols))
    values_whittaker_P[:] = pd.np.nan
    values_whittaker_PR = np.empty((ztime, rows, cols))
    values_whittaker_PR[:] = pd.np.nan
    values_whittaker_PRS = np.empty((ztime, rows, cols))
    values_whittaker_PRS[:] = pd.np.nan
    
    counter = 1
    print ('Running Whittaker-Henderson Smoother...')
    for m,n in product(range(1),range(1)): # rows cols 
        print ('\t{0}/{1}'.format(counter, size_st))
        
        yp = pd.Series(original_P_values[:, m, n])
        yp_inter = yp.interpolate(method='linear')
        yr_p = whfilter(yp_inter, weights, lamb)
        values_whittaker_P[:, m, n] = yr_p

        ypr = pd.Series(original_PR_values[:, m, n])
        ypr_inter = ypr.interpolate(method='linear')
        yr_pr = whfilter(ypr_inter, weights, lamb)
        values_whittaker_PR[:, m, n] = yr_pr
        
        yprs = pd.Series(original_PRS_values[:, m, n])
        yprs_inter = yprs.interpolate(method='linear')
        yr_prs = whfilter(yprs_inter, weights, lamb)
        values_whittaker_PRS[:, m, n] = yr_prs
        
        counter += 1
        
    ncds.variables['whittaker_P_%s_values'%viIndexName][:] = values_whittaker_P*10000
    ncds.variables['whittaker_PR_%s_values'%viIndexName][:] = values_whittaker_PR*10000
    ncds.variables['whittaker_PRS_%s_values'%viIndexName][:] = values_whittaker_PRS*10000
    
    # Close netcdf file
    ncds.close()

def check_fit(ncFile,point):
    """
    Function to check All Smoothing methods results
    using Matplotlib in some pixels
    """
    # Read NetCDF File
    ncds = Dataset(ncFile,'r')
    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    
    # Get x and y
    x = point[0]
    y = point[1]
    
    # Get row and col index closest to x, y in the netcdf file
    y_closest = Y.flat[pd.np.abs(Y - y).argmin()]
    x_closest = X.flat[pd.np.abs(X - x).argmin()]

    row = pd.np.where(Y == y_closest)[0][0]
    col = pd.np.where(X == x_closest)[0][0]
    
    # Read values
    viIndexName = os.path.basename(ncFile).split('_')[0]
    # original
    values_o_p = ncds.variables['original_P_%s_values'%viIndexName][:,row,col]/10000
    values_o_pr = ncds.variables['original_PR_%s_values'%viIndexName][:,row,col]/10000
    values_o_prs = ncds.variables['original_PRS_%s_values'%viIndexName][:,row,col]/10000
    # Hants
    values_h_p = ncds.variables['hants_P_%s_values'%viIndexName][:,row,col]/10000
    values_h_pr = ncds.variables['hants_PR_%s_values'%viIndexName][:,row,col]/10000
    values_h_prs = ncds.variables['hants_PRS_%s_values'%viIndexName][:,row,col]/10000
    combined_h_p = ncds.variables['combined_P_%s_values'%viIndexName][:,row,col]/10000
    combined_h_pr = ncds.variables['combined_PR_%s_values'%viIndexName][:,row,col]/10000
    combined_h_prs = ncds.variables['combined_PRS_%s_values'%viIndexName][:,row,col]/10000
    outliers_h_p = ncds.variables['outliers_P_%s_values'%viIndexName][:,row,col]/10000
    outliers_h_pr = ncds.variables['outliers_PR_%s_values'%viIndexName][:,row,col]/10000
    outliers_h_prs = ncds.variables['outliers_PRS_%s_values'%viIndexName][:,row,col]/10000
    # Savitzky-Golay
#    values_s_p = ncds.variables['savgol_P_%s_values'%viIndexName][:,row,col]/10000
#    values_s_pr = ncds.variables['savgol_PR_%s_values'%viIndexName][:,row,col]/10000
#    values_s_prs = ncds.variables['savgol_PRS_%s_values'%viIndexName][:,row,col]/10000
    # Whittaker
#    values_w_p = ncds.variables['whittaker_P_%s_values'%viIndexName][:,row,col]/10000
#    values_w_pr = ncds.variables['whittaker_PR_%s_values'%viIndexName][:,row,col]/10000
#    values_w_prs = ncds.variables['whittaker_PRS_%s_values'%viIndexName][:,row,col]/10000
    # Time
    times = [pd.to_datetime(Date, format='%Y%m%d') for Date in ncds.variables['time'][:]]
    
    # Figure PlanetScope
    plt.figure() 
    top = 1.15*max(pd.np.nanmax(values_o_p), pd.np.nanmax(values_h_p))#, pd.np.nanmax(combined_h_p), pd.np.nanmax(outliers_h_p))# ,pd.np.nanmax(values_s_p), pd.np.nanmax(values_w_p))
    bottom = min(pd.np.nanmin(values_o_p), pd.np.nanmin(values_h_p))#, pd.np.nanmin(combined_h_p), pd.np.nanmin(outliers_h_p))#   , pd.np.nanmin(values_h_p), pd.np.nanmin(values_w_p)) - 0.05
    ylim = [bottom, top]
        # Plot
    plt.plot(times, values_h_p, 'r-', label='HANTS')
#    plt.plot(times, combined_h_p, 'c--', label='Combined')
#    plt.plot(times, outliers_h_p, 'm.', label='outliers')
#    plt.plot(times, values_s_p, 'g-', label='Savitzky-Golay Filter')
#    plt.plot(times, values_w_p, 'y-', label='Whittaker-Henderson Smoother')
    plt.plot(times, values_o_p, 'b.', label='Original data')
    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: X {0:.2f}, Y {1:.2f} \n PlanetScope Time Serie'.format(x,y))
    plt.axes().set_aspect(0.5*(times[-1] - times[0]).days/(ylim[1] - ylim[0]))
    plt.show()
    
    # Figure PlanetScope & RapidEye
    plt.figure() 
    top = 1.15*max(pd.np.nanmax(values_o_pr), pd.np.nanmax(values_h_pr))#,pd.np.nanmax(combined_h_pr), pd.np.nanmax(outliers_h_pr))# ,pd.np.nanmax(values_s_pr), pd.np.nanmax(values_w_pr))
    bottom = min(pd.np.nanmin(values_o_pr), pd.np.nanmin(values_h_pr))#,pd.np.nanmin(combined_h_pr), pd.np.nanmin(outliers_h_pr))#   , pd.np.nanmin(values_h_pr), pd.np.nanmin(values_w_pr)) - 0.05
    ylim = [bottom, top]
        # Plot
    plt.plot(times, values_h_pr, 'r-', label='HANTS')
#    plt.plot(times, combined_h_pr, 'c--', label='Combined')
#    plt.plot(times, outliers_h_pr, 'm.', label='outliers')
#    plt.plot(times, values_s_pr, 'g-', label='Savitzky-Golay Filter')
#    plt.plot(times, values_w_pr, 'y-', label='Whittaker-Henderson Smoother')
    plt.plot(times, values_o_pr, 'b.', label='Original data')
    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: X {0:.2f}, Y {1:.2f} \n PlanetScope and RapidEye Time Serie'.format(x,y))
    plt.axes().set_aspect(0.5*(times[-1] - times[0]).days/(ylim[1] - ylim[0]))
    plt.show()
    
    # Figure PlanetScope, RapidEye and Sentinel-2
    plt.figure() 
    top = 1.15*max(pd.np.nanmax(values_o_prs), pd.np.nanmax(values_h_prs))#,pd.np.nanmax(combined_h_prs), pd.np.nanmax(outliers_h_prs))# ,pd.np.nanmax(values_s_prs), pd.np.nanmax(values_w_prs))
    bottom = min(pd.np.nanmin(values_o_prs), pd.np.nanmin(values_h_prs))#,pd.np.nanmin(combined_h_prs), pd.np.nanmin(outliers_h_prs))#   , pd.np.nanmin(values_h_prs), pd.np.nanmin(values_w_prs)) - 0.05
    ylim = [bottom, top]
        # Plot
    plt.plot(times, values_h_prs, 'r-', label='HANTS')
#    plt.plot(times, combined_h_prs, 'c--', label='Combined')
#    plt.plot(times, outliers_h_prs, 'm.', label='outliers')
#    plt.plot(times, values_s_prs, 'g-', label='Savitzky-Golay Filter')
#    plt.plot(times, values_w_prs, 'y-', label='Whittaker-Henderson Smoother')
    plt.plot(times, values_o_prs, 'b.', label='Original data')
    plt.ylim(ylim[0], ylim[1])
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('values')
    plt.gcf().autofmt_xdate()
    plt.axes().set_title('Point: X {0:.2f}, Y {1:.2f} \n Entire Time Serie'.format(x,y))
    plt.axes().set_aspect(0.5*(times[-1] - times[0]).days/(ylim[1] - ylim[0]))
    plt.show()

    # Close netcdf file
    ncds.close()
    

if __name__=="__main__":
    
    inPath = "D:/Stage/Output/INDICES/NDVI"
    inVectorFile = "D:/Stage/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"
#
#    ncFile = create_time_series(inPath,inVectorFile)
    ncFile = 'D:/Stage/Output/INDICES/NDVI/NDVI_TIME_SERIES.nc'
#    smooth_hants(ncFile)
#    smooth_savgol(ncFile)
#    smooth_whittaker(ncFile)
    
    Points = [(336142.1, 1602889.5), (335329, 1603060), (334582.4, 1600219.6), (334578.5,1600143.5), (336082.9,1603694.0), \
          (336746.3,1603706.7), (336423.9,1603248.0), (336858.4,1603152.2), (337457.6,1603508.4), (337067.1,1603144.3)]
    
    for point in Points :
        check_fit(ncFile, point)