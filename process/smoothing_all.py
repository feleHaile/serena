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
from scipy.signal import savgol_filter
import scipy as sp
from scipy import sparse, linalg

def create_time_series (inFolder, inVectorFile):
    """
    - Function to create Vegetation Index Time Series in NetCDF4 Format 
    - 3 Time Series variables are created : PlanetScope - PlanetScope & RapidEye - PlanetScope, RapidEye & Sentinel-2
    - Input should be Folder containing GeoTiff files
    - VI Index Values have 10 000 scale factor
    """
    # Dates for Sentinel-2A Images deleting
    S2_Dates = ["20171005","20171015","20171025"]
    
    # Index Folders
    indexFolders = [os.path.join(inFolder,folder) for folder in os.listdir(inFolder) if os.path.isdir(os.path.join(inFolder,folder))]
#    print (indexFolders)
    
    inPath = indexFolders[0]
     
    # Create Date List
    lstFiles = sorted(glob.glob(inPath+'/*.tif'))
    dicFile = {}
    for File in lstFiles :
#        print (File)
        dicFile.update({os.path.basename(File).split('_')[1]:os.path.basename(File).split('_')[2]})
#    print (dicFile)
    origin_time = list(dicFile.keys())
               
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
    outFolder = os.path.join(os.path.dirname(inFolder),'TS')
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    outFile = os.path.join(outFolder,'TIME_SERIES_2.nc')
    ncds = Dataset(outFile ,'w',format='NETCDF4')
#    print (ncds)
    
    # Create Dimensions 
    print ('Creating Dimensions')
    ncds.createDimension("time",len(lstDate))
    ncds.createDimension("X", len(lstX))
    ncds.createDimension("Y", len(lstY))
    ncds.createDimension("time_origin", len(origin_time))
    
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
    
    time_origin_var = ncds.createVariable('time_origin', 'l', ('time_origin'),fill_value=fill_val)
    time_origin_var.standard_name = 'Origin Times in Time Serie'
    time_origin_var.calendar = 'gregorian'
    
    print ('Load General Variables')
    # Load data
    y_var[:] = lstY
    x_var[:] = lstX
    time_var[:] = lstDate
    time_origin_var[:] = origin_time
        
    xOffset = int((ulx- originX) / cellsize)
    yOffset = int((uly - originY) / -cellsize)
    xsize = int((lrx - ulx)/cellsize)  # Raster xsize col
    ysize = int((uly - lry)/cellsize) # Raster ysize lines
    empty_vec = np.empty((len(lstY),len(lstX)))
    empty_vec[:] = -9999.
    
    for j in range (len(indexFolders)) : 
        viIndexName = os.path.basename(indexFolders[j])
        inPath = indexFolders[j]
        lstFiles = sorted(glob.glob(inPath+'/*.tif'))
        FileName = os.path.basename(lstFiles[0]).split('_')[0]+'_%s_%s_'+\
           os.path.basename(lstFiles[0]).split('_')[3]+'_'+os.path.basename(lstFiles[0]).split('_')[4]
#       print (FileName)
    
        vi_var_PRS = ncds.createVariable('original_PRS_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
        vi_var_PRS.long_name = 'PlanetScope, RapidEye & Sentinel-2 %s Index values'%viIndexName
        
        vi_var_cor = ncds.createVariable('original_PRScor_%s_values'%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
        vi_var_cor.long_name = 'PlanetScope, RapidEye & Sentinel-2 %s Index values excluding S2 20171005, 20171015, 20171025 times'%viIndexName
        
        print ('Load %s Variables'%viIndexName)
        
        i = 0
        for Date in lstDate :
            print (Date)
            if (Date in dicFile) : # and dicFile[Date]!='S2')  : # Exclude Sentinel-2 Imagery         
                File = os.path.join(inPath,FileName%(Date,dicFile[Date]))
                print (File)
                with rasterio.open(File,'r') as ds:
                    band = ds.read(1, window=((yOffset, yOffset+ysize),(xOffset, xOffset+xsize)))
                    array = band * 10000
                    array[np.isnan(array)] = fill_val
        #            print (array)
                    vi_var_PRS [i,:,:] = array
                    print ("PRS")
                    if (Date not in S2_Dates):
                        vi_var_cor[i,:,:] =  array
                        print ("PRScor")
                    i+=1
            else :
                vi_var_PRS [i,:,:] =  empty_vec
                vi_var_cor[i,:,:] =  empty_vec
                i+=1
#            
    ncds.close()
    print ('###### NetCDF file created ######')
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
    
def smooth_hants (ncFile, variable, nb=365, nf=3, HiLo='Lo', low=-0.3, high=1, fet=0.05, dod=1, delta=0.25, fill_val=-9999.):
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
    viIndexName = variable.split('_')[2]
    
    # Get Existing Values
    time_var = ncds.variables['time'][:]
    original_values = ncds.variables[variable][:]
    original = original_values/10000
    
    # Create HANTS variables
    hants_varName = 'hants_'+variable.split('_')[1]+'_%s_'+variable.split('_')[3]
    hants_var = ncds.createVariable(hants_varName%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    hants_var.long_name = 'HANTS Smoothing on %s'%ncds.variables[variable].long_name

    [ztime ,rows, cols] = original.shape
    size_st = cols*rows
    values_hants = np.empty((ztime, rows, cols))
    values_hants[:] = pd.np.nan
    
    # Additional parameters
    ni = len(time_var)
    ts = range(ni)

    # Loop
    counter = 1
    print ('Running HANTS...')
    for m,n in product(range(rows),range(cols)):
        print ('\t{0}/{1}'.format(counter, size_st))
    
        y = pd.np.array(original[:, m, n])
        y[pd.np.isnan(y)] = fill_val
        [yr, outliers] = HANTS(ni, nb, nf, y, ts, HiLo, low, high, fet, dod, delta, fill_val)
        values_hants[:, m, n] = yr
#        outliers_hants[:, m, n] = outliers

        counter = counter + 1
    
    values_hants = values_hants * 10000
    ncds.variables[hants_varName%viIndexName][:] = values_hants
    # Close netcdf file
    ncds.close()

# =============================================================================
# Smoothing with Whittaker Smoother
# =============================================================================
def whitsmw (y, w, lamb, d=2):
    """
    Whittaker smoother with weights
    Input:
       y: data series, sampled at equal intervals (arbitrary values allowed when missing, but not NaN!)
       w: weights (0/1 for missing/non-missing data)
       lamb: smoothing parameter; large lambda gives smoother result
       d: order of differences (default = 2)
    Output:
       z: smoothed series
       cve: RMS leave-one-out prediction error
       h: diagonal of hat matrix
    
    Remark: the computation of the hat diagonal for m > 100 is experimental;
    with many missing observation it may fail.
    
    Paul Eilers, 2003
    
    For translation : https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
    http://www.courspython.com/tableaux-numpy.html
    """
    
    # Smoothing
    
    m = sp.size(y)
    E = sparse.eye(m).toarray()
    W = sparse.spdiags(w,0,m,m).toarray()
    D = sp.diff(E, d, axis=0)
    C = linalg.cholesky(W + lamb * (D.conj().T).dot(D))
    z = linalg.solve(C,linalg.solve(C.conj().T,w*y))
    
    return z


def smooth_whittaker(ncFile, variable, lamb=5000, d=2):
    """
    Smooth Time Series in NetCDF4 Format with Whittaker Smoother
    """
    # Read netcdfs
    ncds = Dataset(ncFile, 'r+')
    viIndexName = variable.split('_')[2]
    fill_val = -9999.
    
    # Get Existing Values
    original_values = ncds.variables[variable][:]
    original = original_values/10000
    
    # Create Whittaker Smoother variables
    whit_varName = 'whittaker_'+variable.split('_')[1]+'_%s_'+variable.split('_')[3]
    whit_var = ncds.createVariable(whit_varName%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
    whit_var.long_name = 'Whittaker Smoothing on %s'%ncds.variables[variable].long_name
    
    # Filtering 
    
    [ztime ,rows, cols] = original.shape
    size_st = cols*rows
    values_whit = np.empty((ztime, rows, cols))
    values_whit[:] = pd.np.nan
    
    counter = 1
    print ('Running Whittaker Smoother...')
    for m,n in product(range(rows),range(cols)): # rows cols 
        print ('\t{0}/{1}'.format(counter, size_st))
        
        y = pd.np.array(original[:, m, n])
        ynan = np.where(y==-9999.,np.nan,y)
        w = np.where(y==-9999,0,1) # Negative Values to test : ok
        z = whitsmw (y, w, lamb, d)

        y_new = np.where(y<z,z,y)
#        print (y_new)
        
        di = np.empty((y.shape))
        dmax = 0
        for i in range(y.shape[0]):
            di[i] = np.abs(ynan[i] - z[i])
            if di[i] > dmax:
                dmax = di[i]
#        print (di, dmax)
        
        W = np.where(y<z,1-(di/dmax),1)
#        print (W)
        
        iteration = 1
#        print ('iteration %s'%iteration)
        Fk = Fkmin = 1000000.
        ziter = np.copy(y_new)
        while (Fk <= Fkmin) :
            ziter = whitsmw (ziter, w, lamb=2500)
            lst = []
            for k in range(ziter.shape[0]):
                lst.append(np.abs(ziter[k]-ynan[k])*W[k])
            array = np.array(lst)
            Fk = np.nansum(array)
            if Fk < Fkmin :
                Fkmin = Fk
            else:
                break
            ziter = np.where(y<ziter,ziter,y)
            iteration+=1

        values_whit[:, m, n] = ziter
        counter += 1

    values_whit = values_whit * 10000
    ncds.variables[whit_varName%viIndexName][:] = values_whit #  à remettre

    # Close netcdf file
    ncds.close()

# =============================================================================
# Smoothing with Savitzky-Golay Filter
# =============================================================================

def smooth_savgol(ncFile, variable, Points, lstWindow_size, lstOrder):
    """
    Read Time Series variable in NetCDF4 Format and Smooth values with Savitzky-Golay 
    Savitzky-Golay parameters are :
       windows_length : Size of Filter Window (must be odd 2m+1 with m as window half-width 
       polyorder : Order of Fit Polynomial Function
    Save Smoothed Values for current Point and setted parameters in CSV File
    """
    
    # Read netcdfs
    ncds = Dataset(ncFile, 'r')
    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    times = ncds.variables['time'][:]
    
    viIndexName = variable.split('_')[2]
    outFolder = os.path.join(os.path.dirname(ncFile),str(viIndexName)+'_aggregate')
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    outCSV = os.path.join(outFolder,'pxl_%s_savgol_%s.csv'%(viIndexName,variable.split('_')[1]))
    
    csvDict = {}
    
    for point in Points:
        # Get Row and Col number for Point
        id_point = point[0]
        x_point = point[1]
        y_point = point[2]
        
        # Get row and col index closest to x, y in the netcdf file
        y_closest = Y.flat[pd.np.abs(Y - y_point).argmin()]
        x_closest = X.flat[pd.np.abs(X - x_point).argmin()]
    
        row = pd.np.where(Y == y_closest)[0][0]
        col = pd.np.where(X == x_closest)[0][0]
    
        # Get Existing Values
        original_values = ncds.variables[variable][:,row,col]
        original = original_values/10000
        
        for order in lstOrder : 
            for window_size in lstWindow_size:
                # Filtering
                y = pd.np.array(original) # m, n
                yn = np.where(y==-9999.,np.nan,y)
                y_inter = pd.Series(yn).interpolate(method='linear')
                ys = savgol_filter(y_inter.values,window_size, order)
                
                for i in range(len(times)):
                    
                    csvDict.setdefault("Date",[]).append(times[i])
                    csvDict.setdefault("Original",[]).append(yn[i])
                    csvDict.setdefault("Interpolate",[]).append(y_inter[i])
                    csvDict.setdefault("Savgol",[]).append(ys[i])
                
                
                    csvDict.setdefault("window",[]).append(window_size)
                    csvDict.setdefault("order",[]).append(order)
                    csvDict.setdefault("Point_id",[]).append(id_point)
                    csvDict.setdefault("Point_X",[]).append(x_point)
                    csvDict.setdefault("Point_Y",[]).append(y_point)
                    
#    print (csvDict)
    outdf = pd.DataFrame.from_dict(csvDict)
    outdf.to_csv(outCSV,index=False)

    # Close netcdf file
    ncds.close()


if __name__=="__main__":
    
# =============================================================================
#     NDVI
# =============================================================================
    
#    inFolder = "E:/Stage2018/Output/INDICES"
#    inVectorFile = "E:/Stage2018/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"

#    ncFile = create_time_series(inFolder,inVectorFile)
    ncFile = 'D:/TIME_SERIES.nc'
#    
    lstVar = ["original_PRScor_NDVI_values","original_PRS_NDVI_values"]
#              "original_PRScor_MSAVI2_values", "original_PRS_MSAVI2_values"]
#    
#    for variable in lstVar:
#        print (variable)
##        smooth_hants(ncFile, variable)
#        smooth_whittaker(ncFile, variable)
    
    Points = [(1,336142.1, 1602889.5), (2,335329, 1603060), (3,334582.4, 1600219.6), (4,334578.5,1600143.5), (5,336082.9,1603694.0), \
          (6,336746.3,1603706.7), (7,336423.9,1603248.0), (8,336858.4,1603152.2), (9,337457.6,1603508.4), (10,337067.1,1603144.3)]

    lstOrder = [2,3,4]
    lstWindow_size = [9,11,13,15]
    
    for variable in lstVar:
    #    print(variable,point,order,window_size)
        smooth_savgol(ncFile, variable, Points, lstWindow_size, lstOrder)
                    
#    ncFile = "E:/Stage2018/Output/TS/TIME_SERIES_2.nc"
#    variable = "original_PRS_NDVI_values"
#    smooth_savgol(ncFile, variable, point=Points[0], window_size=11, order=4)

