# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 16:21:59 2018

@author: User
"""
import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from itertools import product
from copy import deepcopy
import math
import scipy as sp
from scipy import sparse, linalg

# =====================
# Smoothing with HANTS
# =====================
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

# ==================================
# Smoothing with Whittaker Smoother
# ==================================
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
    
def smoothing (inCSV, method):
    """
    Smooth Time Series in CSV with Harmonic Analysis of Time series or Whittaker Smoother
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
    data = pd.read_csv(inCSV)
    # print (data)
    
    viIndexName = os.path.basename(inCSV).split('_')[1]
    # print (viIndexName)
    
    times = pd.date_range(start='2017-05-08', end = '2017-11-19', freq='D')
    
    outDict = {}
    
    if method == "HANTS":
        nb=365
        nf=3
        HiLo='Lo'
        low=-0.3
        high=1
        fet=0.05
        dod=1
        delta=0.25
        fill_val=-9999
        
        # Additional parameters
        ni = len(times)
        ts = range(ni)

    for i in range(len(data)): # Each Plot
        print ("ID: ",data.ID[i])
        # Append each date value to list
        lstValues = []
        for j in range(len(times)): # dates
            date = str(times.strftime('%Y%m%d')[j])+'Mean'
            if date in data.columns : 
                lstValues.append(data[date][i])
            else :
                lstValues.append(-9999.)
        
        # Create Dataframe to extract metrics for current plot
        plot_df = pd.DataFrame({"Date":times,"Value":lstValues})
        # print (plot_df)
        y = pd.np.array(plot_df.Value)
        
        if method == "HANTS":
            y[pd.np.isnan(y)] = fill_val
            [yr, outliers] = HANTS(ni, nb, nf, y, ts, HiLo, low, high, fet, dod, delta, fill_val)
            
            outDict.setdefault("ID",[]).append(data.ID[i])
            
            for j in range(len(times)):
                date = str(times.strftime('%Y%m%d')[j])+'Mean'
                outDict.setdefault(date,[]).append(yr[j])
    
    outCSV = os.path.join("H:/Stage2018/Process/LISSAGE/GDVI_aggregate_notree/agg_GDVI_hants_PRScor.csv")
    outdf = pd.DataFrame.from_dict(outDict)
    outdf.to_csv(outCSV,index=False)
            

    
#    
#    
#    # Get Existing Values
#    time_var = ncds.variables['time'][:]
#    original_values = ncds.variables[variable][:]
#    original = original_values/10000
#    
#    # Create HANTS variables
#    hants_varName = 'hants_'+variable.split('_')[1]+'_%s_'+variable.split('_')[3]
#    hants_var = ncds.createVariable(hants_varName%viIndexName, 'i2', ('time', 'Y', 'X'), fill_value=fill_val)
#    hants_var.long_name = 'HANTS Smoothing on %s'%ncds.variables[variable].long_name
#
#    [ztime ,rows, cols] = original.shape
#    size_st = cols*rows
#    values_hants = np.empty((ztime, rows, cols))
#    values_hants[:] = pd.np.nan
#    
#    
#
#    # Loop
#    counter = 1
#    print ('Running HANTS...')
#    for m,n in product(range(rows),range(cols)): # rows cols
#        print ('\t{0}/{1}'.format(counter, size_st))
#    
#        y = pd.np.array(original[:, m, n]) #m n 
#        y[pd.np.isnan(y)] = fill_val
#        
#        values_hants[:, m, n] = yr
##        outliers_hants[:, m, n] = outliers
#
#        counter = counter + 1
#    
#    values_hants = values_hants * 10000
#    ncds.variables[hants_varName%viIndexName][:] = values_hants
#    # Close netcdf file
#    ncds.close()


if __name__=="__main__":
    
    inCSV = "H:/Stage2018/Process/LISSAGE/GDVI_aggregate_notree/agg_GDVI_original_PRScor.csv"
    smoothing(inCSV, method="HANTS")
    