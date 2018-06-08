#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 12:44:56 2018

@author: je

Aggregate VI Index Pixels Values in Parcels Mean Value (Ground Database Vector File) 
to check Smoothing Results

"""
#import rasterio, fiona
import gdal, ogr, osr
from osgeo.gdalconst import GA_ReadOnly
#import geopandas as gpd
#from shapely.wkt import loads
import pandas as pd
import os
import numpy as np
from netCDF4 import Dataset
import sys
#from pprint import pprint

def create_csv(outCSV, time, suffix, xsize, ysize, originX, originY, Proj, rGeoT, v_values) :
    """
    Aggregate Parcels Mean Values and Write them in CSV
    
    """
    # Read ShapeFile 
    shpds = ogr.Open(inVectorFile,GA_ReadOnly)
    if shpds is None:
        print ('Couldn\'t open ShapeFile. Please Check it.')
        sys.exit (1)
#        print (shpds)
    layer = shpds.GetLayer()
    
    # Create in Memory Gdal and Ogr Driver 
    rmem_drv = gdal.GetDriverByName("MEM")
    vmem_drv = ogr.GetDriverByName("MEMORY")
       
    # For each Feature, Get Geometry, attributes and Compute Zonal Stats
    csvDict = {}
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        print ("Feature %s ID %s"%(i,feature.GetField('ID')))
        geometry = feature.GetGeometryRef()
        extent = geometry.GetEnvelope()
        
        csvDict.setdefault('ID',[]).append(feature.GetField('ID'))
        
#        print (extent)
        ulx = extent[0]
        uly = extent[3]
        lrx = extent[1]
        lry = extent[2]
        
        vGeoT = (ulx,3.,0,uly,0,-3.)
        cellsize = 3.
        
        xOffset = int((ulx - originX) / cellsize)
        yOffset = int((uly - originY) / -cellsize)

        col_nb = int((lrx - ulx)/cellsize) # Raster xsize col
        row_nb = int((uly - lry)/cellsize)  # Raster ysize lines
        
#        lstMean = []
        counter = 0
        for Date in time : 
#            print (Date)
            # Create in memory Raster for current date
            dest = rmem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
            dest.SetGeoTransform(rGeoT)
            dest.SetProjection(Proj.ExportToWkt())
            dest.GetRasterBand(1).WriteArray(v_values[counter,:,:])

            valuesArray = dest.ReadAsArray(xOffset,yOffset,col_nb,row_nb)
            
            # Create Raster in Memory for feature Rasterization
            rfeat_ds = rmem_drv.Create('', col_nb, row_nb, 1, gdal.GDT_Byte)
            rfeat_ds.SetGeoTransform(vGeoT)
            rfeat_ds.SetProjection(Proj.ExportToWkt())
            
            # Create a temporary vector layer in memory containing current feature
            vfeat_ds = vmem_drv.CreateDataSource('Out')
            vfeat_layer = vfeat_ds.CreateLayer('Polygon', Proj, ogr.wkbPolygon)
            vfeat_layer.CreateFeature(feature.Clone())
            # Rasterize feature
            gdal.RasterizeLayer(rfeat_ds, [1], vfeat_layer, burn_values=[1])
            
            # Now check raster feature pixel value to select valid extent pixels for Zonal Stats
                # Mask the source data array with our current feature
                # we take the logical_not to flip 0<->1 to get the correct mask effect
                # we also mask out nodata values explictly
            rfeatArray = rfeat_ds.ReadAsArray()
            nodata = np.nan
            masked = np.ma.MaskedArray(valuesArray,mask=np.logical_or(valuesArray == nodata,np.logical_not(rfeatArray)))
#            print (masked)
            
            mean = masked.mean()
            std = masked.std()
#            masked.max()
#            masked.sum()
#            masked.count()
            csvDict.setdefault(str(Date)+suffix%'Mean',[]).append(mean)
            csvDict.setdefault(str(Date)+suffix%'Std',[]).append(std)
    
            rfeat_ds = None
            vfeat_ds = None
            
            counter += 1
    
    # Create CSV with pandas from dict
    outdf = pd.DataFrame.from_dict(csvDict)
    outdf.to_csv(outCSV,index=False)

    rmem_drv = None
    vmem_drv = None
    shpds = None
    
def aggregate (inFile, inVectorFile, variable) :
    """
    Agregate inFile Raster Pixels Values in Parcel (Vector File) Mean Value
    Based on perrygeo/zonal_stats.py : https://gist.github.com/perrygeo/5667173
    and gdal/ogr cookbook : https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#calculate-zonal-statistics
    Raster and Vector Files must have Same SRS
    
    """
    # Read NetCDF File 
    # Extract some parameters for Zonal Stats
    ncds = Dataset(inFile)
    viIndexName = os.path.basename(inFile).split('_')[0]
    
    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    xsize = len(X)
    ysize = len(Y)
    time = [Date for Date in ncds.variables['time'][:]]
#    print (time)
    Proj = osr.SpatialReference()
    Proj.ImportFromEPSG(32628)
    
    originX = min(X)-0.5*3.
    originY = max(Y)-0.5*3.
    
    rGeoT = (originX,3.,0,originY,0,-3.)
    
#    lstVariables = [variable%viIndexName], 'original_PR_%s_values'%viIndexName, 'original_PRS_%s_values'%viIndexName,
#                 'hants_P_%s_values'%viIndexName, 'hants_PR_%s_values'%viIndexName, 'hants_PRS_%s_values'%viIndexName]
#                 ,'savgol_P_%s_values'%viIndexName, 'savgol_PR_%s_values'%viIndexName, 'savgol_PRS_%s_values'%viIndexName,
#                 'whittaker_P_%s_values'%viIndexName, 'whittaker_PR_%s_values'%viIndexName, 'whittaker_PRS_%s_values'%viIndexName]
    # CSV Header
#    for variable in lstVariables :
    if variable.startswith('original'):
        suffix = '_O%s'
#        time2 = [str(Date)+'_OMean' for Date in time]
#        time2.insert(0,"ID")
    else :
        suffix = '_S%s'
#        time2 = [str(Date)+'_SMean' for Date in time]
#        time2.insert(0,"ID")
#    print (time2)
     
    # Extract Original and Hants Values
    variable_values = np.array(ncds.variables[variable%viIndexName][:])
    v_values = np.where(variable_values==-9999.,np.nan,variable_values)
    v_values = v_values / 10000
    # CSV Name
    outCSV = os.path.join(os.path.dirname(inFile),'Mean_%s_%s_%s_iter.csv'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    create_csv(outCSV, time, suffix, xsize, ysize, originX, originY, Proj, rGeoT, v_values)
    
    # CSVt File for Attribute Join in GIS Software
    outCSVt = os.path.join(os.path.dirname(inFile),'Mean_%s_%s_%s_iter.csvt'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    with open(outCSVt,'w') as csvt_file:
        csvt_file.write('"String",') 
        for i in range(len(time)-1) :
            csvt_file.write('"Real",')
        csvt_file.write('"Real"')
        
    # Close netCDF dataset
    ncds = None

def mean_profile (hantsCSV, whitCSV) :
    """
    Compute mean value about Hants and Whittaker Smoother for each plot
    """
    outFile = os.path.basename(hantsCSV).split('_')[0]+'_'+os.path.basename(hantsCSV).split('_')[1]+'_SmoothAvg_'+os.path.basename(whitCSV).split('_')[3]+'_'+os.path.basename(whitCSV).split('_')[4]
    hants_df = pd.read_csv(hantsCSV)
    whit_df = pd.read_csv(whitCSV)
#    print (hants_values, whit_values)
    df_concat = pd.concat((hants_df, whit_df))
    by_row_index = df_concat.groupby(df_concat.index)
    df_means = by_row_index.mean()
    df_means['ID'] =  hants_df['ID'] # Create ID Column
    df_means.to_csv(os.path.join(os.path.dirname(hantsCSV),outFile),index = False)
    
    # CSVt File for Attribute Join in GIS Software
    outCSVt = os.path.join(os.path.dirname(hantsCSV),outFile+'t')
    with open(outCSVt,'w') as csvt_file:
        csvt_file.write('"String",') 
        for i in range(len(hants_df.columns)-2) : # ID column are not considered
            csvt_file.write('"Real",')
        csvt_file.write('"Real"')

if __name__=="__main__":
    
# =============================================================================
#     NDVI
# =============================================================================
    
    inFile = "D:/Stage/Output/INDICES/NDVI/NDVI_TIME_SERIES.nc"
    inVectorFile = "D:/Stage/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"

#    variables = ["original_PRScor_%s_values","hants_PRScor_%s_values"]#,"whittaker_PRScor_%s_values"]#,"savgol_PRScor_%s_values"]
    
#    variable = "original_PRS_%s_values"
    variable = "whittaker_PRScor_%s_values"
#    aggregate (inFile, inVectorFile, variable)
    
#    for variable in variables :
#        print (variable%os.path.basename(inFile).split('_')[0])
#        aggregate (inFile, inVectorFile, variable)
        
# =============================================================================
#     MSAVI2
# =============================================================================
    
#    inFile = "D:/Stage/Output/INDICES/NDVI/SAVI_TIME_SERIES.nc"
#    inVectorFile = "D:/Stage/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"
#
#    variables = ["original_PRScor_%s_values","hants_PRScor_%s_values","whittaker_PRScor_%s_values","savgol_PRScor_%s_values"]
#    for variable in variables :
#        print (variable%os.path.basename(inFile).split('_')[0])
#        aggregate (inFile, inVectorFile, variable)
    
# =============================================================================
#     Mean Profile
# =============================================================================
    hantsCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_hants_PRScor.csv"
    whitCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_whittaker_PRScor_iter.csv"
    mean_profile (hantsCSV, whitCSV)
