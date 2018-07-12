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
        SystCult = str(join_df.CropSyst[i]+' '+join_df.Crop_1[i])
#        print (SystCult)
        lstValues= [] # List to append values
        
        # Append each date value to list
        for j in range(len(times)):
            date = str(times.strftime('%Y%m%d')[j])+'Mean'
            lstValues.append(join_df[date][i])
        
        # Create Dataframe to extract metrics for current plot
        plot_df = pd.DataFrame({"Date":times,"Value":lstValues})
#        print (plot_df["Date"][1].strftime('%Y%m%d'))
        
    # =============================================================================
    #         First Method : NDVI Ratio (NDVI - NDVImin) / (NDVImax - NDVImin)
    # =============================================================================
        
        ndvi_min = plot_df.Value.min()
        ndvi_max = plot_df.Value.max()
#        print (ndvi_min, ndvi_max)
        
        # Get SOS (10% for Pure and Mixed Millet Plots, 20% for Mixed Groundnut)
        # Get EOS (80% for Pure and Mixed Millet and 50% for Mixed Groundnut)
        
        if (SystCult == "Pure Millet" or SystCult == "Mixed Millet") : 
            
            for k in range(len(times)):
                ratio_sos = (plot_df.Value[k] - ndvi_min)/(ndvi_max-ndvi_min)
                if (ratio_sos >= 0.1):
                    break
            SOS = int(plot_df["Date"][k].strftime('%j'))
            SOS_index = k
#            plot_df["Date"][k].strftime('%Y%m%d')
#            print (ratio, plot_df["Date"][k].strftime('%Y%m%d'))
            
            for k in range(len(times)):
                ratio_eos = (plot_df.Value[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
                if (ratio_eos >= 0.8):
                    break
            EOS = int(plot_df["Date"][len(times)-(k+1)].strftime('%j'))
            EOS_index = len(times)-(k+1)
            
            
        elif SystCult == "Mixed Groundnut" :
            for k in range(len(times)):
                ratio_sos = (plot_df.Value[k] - ndvi_min)/(ndvi_max-ndvi_min)
                if (ratio_sos >= 0.2):
                    break
            SOS = int(plot_df["Date"][k].strftime('%j'))
            SOS_index = k
            
            for k in range(len(times)):
                ratio_eos = (plot_df.Value[len(times)-(k+1)] - ndvi_min)/(ndvi_max-ndvi_min)
                if (ratio_eos >= 0.5):
                    break
            EOS = int(plot_df["Date"][len(times)-(k+1)].strftime('%j'))
            EOS_index = len(times)-(k+1)
        
        # Get LOS
        LOS = EOS - SOS
        
        # Get AMPL
        AMPL = ndvi_max - ndvi_min
        
        # Integrated values from SOS to EOS
        CUM_SOS_EOS = 0
        for m in range(SOS_index,EOS_index+1) :
            CUM_SOS_EOS += plot_df["Value"][m]
            
        # Integrated values from SOS to MAX
        MAX_index = plot_df.loc[plot_df["Value"]==ndvi_max].index[0]
        CUM_SOS_MAX = 0
        for m in range(SOS_index,MAX_index+1) :
            CUM_SOS_MAX += plot_df["Value"][m]
        
        # RATE_SOS_MAX Value
        RATE_SOS_MAX = 100 * (ndvi_max - plot_df["Value"][SOS_index]) / plot_df["Value"][SOS_index]
        
        # RATE_MAX_EOS Value
        RATE_MAX_EOS = 100 * (plot_df["Value"][EOS_index] - ndvi_max) / ndvi_max
        
#        print (SOS, EOS, LOS, ndvi_max, AMPL, CUM_SOS_EOS, RATE_SOS_MAX, RATE_MAX_EOS)
   
        # Add to outDict
        outDict.setdefault('ID',[]).append(join_df.ID[i])
        # Integrated Values from SOS to 100 days after SOS for every 5 days 
        for m in range (20): # 0 à 100/5 avec un pas de 5
            CUM = 0
            for n in range (5*m,5*m+5):
                CUM += plot_df["Value"][SOS_index+n]
            outDict.setdefault('CUM_%s_%s'%(str(5*m),str(5*m+5)),[]).append(CUM)
        
        outDict.setdefault('SOS',[]).append(SOS)
        outDict.setdefault('EOS',[]).append(SOS)
        outDict.setdefault('LOS',[]).append(LOS)
        outDict.setdefault('MAX',[]).append(ndvi_max)
        outDict.setdefault('AMPL',[]).append(AMPL)
        outDict.setdefault('CUM_SOS_EOS',[]).append(CUM_SOS_EOS)
        outDict.setdefault('CUM_SOS_MAX',[]).append(CUM_SOS_MAX)
        outDict.setdefault('RATE_SOS_MAX',[]).append(RATE_SOS_MAX)
        outDict.setdefault('RATE_MAX_EOS',[]).append(RATE_MAX_EOS)
        
#        outDict.setdefault('FieldType',[]).append(join_df.FieldType[i])
#        outDict.setdefault('Projet',[]).append(join_df.Projet[i])
#        outDict.setdefault('SystCult',[]).append(SystCult)
#        if not (join_df.Crop_2[i] is None or join_df.Crop_3[i] is None):
#            outDict.setdefault('MixedCrops',[]).append(str(join_df.Crop_2[i]+' '+join_df.Crop_3[i]))
#        else : 
#            outDict.setdefault('MixedCrops',[]).append(join_df.Crop_2[i])
#
#        outDict.setdefault('Nb_Tree',[]).append(join_df.Projet[i])
#        outDict.setdefault('Yield',[]).append(float(join_df.Rdt[i]))
#        outDict.setdefault('Biomass',[]).append(float(join_df.Biom[i]))
   
    outName = 'METRICS_'+os.path.basename(inCSV).split('_',1)[1]
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
    rows, cols = variable_values.shape[1],variable_values.shape[2]
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
    
    inVector = "/home/je/Bureau/Stage/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor.geojson"
    inCSV = "/home/je/Bureau/Stage/Output/TS/NDVI_aggregate/agg_NDVI_whittaker_PRScor.csv"
    outPath = "/home/je/Bureau/Stage/Output/METRICS"

    extract_metrics_plots (inVector,inCSV,outPath)
    
    # lstCSV = ["D:/Stage/TS/NDVI_aggregate/agg_NDVI_hants_PRScor.csv","D:/Stage/TS/NDVI_aggregate/agg_NDVI_whittaker_PRScor.csv",
    #           "D:/Stage/TS/NDVI_aggregate/agg_NDVI_hants_PRS.csv","D:/Stage/TS/NDVI_aggregate/agg_NDVI_whittaker_PRS.csv"]
    
#    for inCSV in lstCSV :
#        extract_metrics_plots (inVector,inCSV,outPath)

#     ncFile = "D:/Stage/TS/TIME_SERIES.nc"
#     lstVariables = ["hants_PRScor_NDVI_values","whittaker_PRScor_NDVI_values"]
#     for variable in lstVariables :
#        extract_metrics_pixels(ncFile,variable,outPath)
