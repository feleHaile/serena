# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 13:50:00 2018

@author: je

Extract Phenology Metrics Sart of Season (SOS), End of Season (EOS) and others & Estimate Sowing dates

"""

import pandas as pd
import geopandas as gpd
import os
from netCDF4 import Dataset
import numpy as np
from itertools import product
import rasterio 

def extract_metrics_plots (inVector, inCSV, outPath):
    """
    Function to extract phenology metrics SOS, EOS ... using treshlod on NDVI Ratio at Plot Scale
    INPUTS
     - Vector File : GeoJSON 
     - CSV File : NDVI Smoothing Profile (Hants, Whittaker for PRSCor or PRS)
     - outPath : Path to write Phenology metrics CSV Files
    """
    # Read GeoJSOn and CSV Files
    gdf = gpd.read_file(inVector)
    csv_df = pd.read_csv(inCSV)
    
    # Join by attribute ID
    join_df = gdf.merge(csv_df,on='ID',how='left')
#    print (join_df)
    
    # Create List of Time Series dates to get values
    times = pd.date_range(start='2017-05-08', end = '2017-11-19', freq='D')
#    print (times)
    
    # For each plot, get Mean values from 20170508 to 20171119
    outDict = {}
    
    for i in range(len(join_df['ID'])): 
#        print (join_df.ID[i]) # ID
        
        lstValues= [] # List to append values
        
        # Append each date value to list
        for j in range(len(times)):
            date = str(times.strftime('%Y%m%d')[j])+'Mean'
            lstValues.append(join_df[date][i])
        
        # Create Dataframe to extract metrics for current plot
        plot_df = pd.DataFrame({"Date":times,"Values":lstValues})
#        print (plot_df["Date"][1].strftime('%Y%m%d'))
        
# =============================================================================
#         First Method : NDVI Ratio (NDVI - NDVImin) / (NDVImax - NDVImin)
# =============================================================================
        
        ndvi_min = plot_df.Values.min()
        ndvi_max = plot_df.Values.max()
#        print (ndvi_min, ndvi_max)
        
        # Get SOS
        # 10%
        for k in range(len(times)):
            ratio_sos = (plot_df.Values[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.1):
                break
        SOS_10 = plot_df["Date"][k].strftime('%Y%m%d')
#        print (ratio, plot_df["Date"][k].strftime('%Y%m%d'))
        
        # 20%
        for k in range(len(times)):
            ratio_sos = (plot_df.Values[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.2):
                break
        SOS_20 = plot_df["Date"][k].strftime('%Y%m%d')
        
        # 30%
        for k in range(len(times)):
            ratio_sos = (plot_df.Values[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.3):
                break
        SOS_30 = plot_df["Date"][k].strftime('%Y%m%d')
            
        # 50%
        for k in range(len(times)):
            ratio_sos = (plot_df.Values[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.5):
                break
        SOS_50 = plot_df["Date"][k].strftime('%Y%m%d')
        
        
        # Get EOS
        # 50%
        for k in range(len(times)):
            ratio_eos = (plot_df.Values[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.5):
                break
#        print (ratio, plot_df["Date"][k].strftime('%Y%m%d'))
        EOS_50 = plot_df["Date"][len(times)-(k+1)].strftime('%Y%m%d')
        
        # 60%
        for k in range(len(times)):
            ratio_eos = (plot_df.Values[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.6):
                break
        EOS_60 = plot_df["Date"][len(times)-(k+1)].strftime('%Y%m%d')
        
        # 70%
        for k in range(len(times)):
            ratio_eos = (plot_df.Values[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.7):
                break
        EOS_70 = plot_df["Date"][len(times)-(k+1)].strftime('%Y%m%d')
        
        # 80%
        for k in range(len(times)):
            ratio_eos = (plot_df.Values[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.8):
                break
        EOS_80 = plot_df["Date"][len(times)-(k+1)].strftime('%Y%m%d')
        
#        print (join_df.ID[i], SOS, join_df.DateSemi[i], EOS, join_df.DateReco[i])
        
        # Add to outDict
        outDict.setdefault('ID',[]).append(join_df.ID[i])
        outDict.setdefault('SOS_10',[]).append(SOS_10)
        outDict.setdefault('SOS_20',[]).append(SOS_20)
        outDict.setdefault('SOS_30',[]).append(SOS_30)
        outDict.setdefault('SOS_50',[]).append(SOS_50)
        
        outDict.setdefault('EOS_50',[]).append(EOS_50)
        outDict.setdefault('EOS_60',[]).append(EOS_60)
        outDict.setdefault('EOS_70',[]).append(EOS_70)
        outDict.setdefault('EOS_80',[]).append(EOS_80)

        outDict.setdefault('Semi',[]).append(join_df.DateSemi[i])
        outDict.setdefault('Recolte',[]).append(join_df.DateReco[i])
        outDict.setdefault('Projet',[]).append(join_df.Projet[i])
        outDict.setdefault('SystCult',[]).append(join_df.CropSyst[i]+'_'+join_df.Crop_1[i])
    
    outName = 'metrics_'+os.path.basename(inCSV).split('_',1)[1]
    outCSV = os.path.join(outPath,outName)
    outdf = pd.DataFrame.from_dict(outDict)
    outdf.to_csv(outCSV,index=False)
    

def extract_metrics_pixels(ncFile,variable,outPath):
    """
    Function to extract phenology metrics SOS, EOS ... using treshlod on NDVI Ratio at pixel scale
    INPUTS
     - netCDF4 File containing smoothing values
     - variable : hants/whittaker_PRScor/PRS values to consider
     - outPath : Path to write phenology metrics Tiff Files
    """
    
    # Read and Get netCDF4 File values

    ncds = Dataset(ncFile)

    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    times = [pd.to_datetime(Date, format='%Y%m%d') for Date in ncds.variables['time'][:]]
    
    variable_values = ncds.variables[variable][:]
    variable_values = variable_values/10000
#    print (time)
    [ztime ,rows, cols] = variable_values.shape
    size_st = cols*rows
    
    SOS_values = np.empty((4, rows, cols),dtype=np.uint16) # 10,20,30,50
    
    EOS_values = np.empty((4, rows, cols),dtype=np.uint16) # 50,60,70,80

    
    # For each pixel, get SOS and EOS in Julian date value using different tresholds 
    counter = 1
    print ('Extracting Phenology Metrics...')
    for m,n in product(range(rows),range(cols)): # rows cols 
        print ('\t{0}/{1}'.format(counter, size_st))
        
        y = pd.np.array(variable_values[:, m, n]) # m,n
        
        ndvi_min = min(y)
        ndvi_max = max(y)
#        print (ndvi_min,ndvi_max)
        
        # SOS 
        # 10%
        for k in range(len(times)):
            ratio_sos = (y[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.1):
                break
        SOS_values[0,m,n] = int(times[k].strftime('%j'))
#        print (ratio, plot_df["Date"][k].strftime('%Y%m%d'))
        
        # 20%
        for k in range(len(times)):
            ratio_sos = (y[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.2):
                break
        SOS_values[1,m,n] = int(times[k].strftime('%j'))
        
        # 30%
        for k in range(len(times)):
            ratio_sos = (y[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.3):
                break
        SOS_values[2,m,n] = int(times[k].strftime('%j'))
            
        # 50%
        for k in range(len(times)):
            ratio_sos = (y[k] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_sos >= 0.5):
                break
        SOS_values[3,m,n] = int(times[k].strftime('%j'))
        
        
        # EOS
        # 50%
        for k in range(len(times)):
            ratio_eos = (y[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.5):
                break
#        print (ratio, plot_df["Date"][k].strftime('%Y%m%d'))
        EOS_values[0,m,n] = int(times[len(times)-(k+1)].strftime('%j'))
        
        # 60%
        for k in range(len(times)):
            ratio_eos = (y[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.6):
                break
        EOS_values[1,m,n] = int(times[len(times)-(k+1)].strftime('%j'))
        
        # 70%
        for k in range(len(times)):
            ratio_eos = (y[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.7):
                break
        EOS_values[2,m,n] = int(times[len(times)-(k+1)].strftime('%j'))
        
        # 80%
        for k in range(len(times)):
            ratio_eos = (y[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
#            print (ratio)
            if (ratio_eos >= 0.8):
                break
        EOS_values[3,m,n] = int(times[len(times)-(k+1)].strftime('%j'))
    
        counter+=1
    
    profile = {'driver': 'GTiff', 'dtype': 'uint16', 'nodata': 0, 'width': len(X), 'height': len(Y), 'count': 4, 
               'crs': ({'init': 'epsg:32628'}), 'transform': (min(X), 3.0, 0.0, max(Y), 0.0, -3.0), 
               'affine': rasterio.Affine(3.0, 0.0, min(X),0.0, -3.0, max(Y)), 'tiled': False, 'compress': 'lzw', 'interleave': 'band'}
        
    outName = variable.split('_')[2] + '_%s_' + variable.split('_')[0]+'_'+variable.split('_')[1]+'.tif'
    
    outFile_sos = os.path.join(outPath,outName%"SOS")
    outFile_eos = os.path.join(outPath,outName%"EOS")
    
    with rasterio.open(outFile_sos,'w', **profile) as outDS :
        for i in range(4):
            outDS.write(SOS_values[i,:,:],i+1)   
    
    with rasterio.open(outFile_eos,'w', **profile) as outDS :
        for i in range(4):
            outDS.write(EOS_values[i,:,:],i+1)
  
if  __name__=='__main__':
    
    inVector = "D:/Stage/TS/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor.geojson"
    inCSV = "D:/Stage/TS/NDVI_aggregate/agg_NDVI_whittaker_PRScor.csv"
    
    lstCSV = ["D:/Stage/TS/NDVI_aggregate/agg_NDVI_hants_PRScor.csv","D:/Stage/TS/NDVI_aggregate/agg_NDVI_whittaker_PRScor.csv",
              "D:/Stage/TS/NDVI_aggregate/agg_NDVI_hants_PRS.csv","D:/Stage/TS/NDVI_aggregate/agg_NDVI_whittaker_PRS.csv"]
    
    outPath = "D:/Stage/Metrics"
    
#    for inCSV in lstCSV :
#        extract_metrics_plots (inVector,inCSV,outPath)

    ncFile = "D:/Stage/TS/TIME_SERIES.nc"
    lstVariables = ["hants_PRScor_NDVI_values","whittaker_PRScor_NDVI_values"]
#    for variable in lstVariables :
#        extract_metrics_pixels(ncFile,variable,outPath)
