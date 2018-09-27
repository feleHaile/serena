# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 17:06:31 2018

@author: JE
"""

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
from itertools import product
#from pprint import pprint

def create_csv(outCSV, time, origin_time, xsize, ysize, originX, originY, Proj, rGeoT, v_values, ovalues_20170618) :
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
        
        # Mask Tree Value
        # Create in memory Raster for date 18-06-2017
        
        dest = rmem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
        dest.SetGeoTransform(rGeoT)
        dest.SetProjection(Proj.ExportToWkt())
        dest.GetRasterBand(1).WriteArray(ovalues_20170618)
        dest.GetRasterBand(1).SetNoDataValue(np.nan)
        
        plot_array = dest.GetRasterBand(1).ReadAsArray(xOffset,yOffset,col_nb,row_nb)
        # print (plot_array)

        tree_mask_array = np.ones((row_nb,col_nb),dtype=int)
        tree_mask_array = np.where(plot_array>0.16,0,plot_array)
        dest = None
        # print (tree_mask_array)
#        driver = gdal.GetDriverByName("GTiff")
#        ds = driver.Create("/home/je/Bureau/Stage/Output/MASK_TREE/Feature_%s_ID_%s.tif"%(i,feature.GetField('ID')), col_nb, row_nb, 1, gdal.GDT_Byte) # 7 is Float_64 Gdal datatype
#        ds.SetGeoTransform(vGeoT)
#        ds.SetProjection(Proj.ExportToWkt())
#        ds.GetRasterBand(1).WriteArray(tree_mask_array)
        
        counter = 0
        for Date in time : 
            if (os.path.basename(outCSV).split('_')[2]=='original' and Date in origin_time) or (os.path.basename(outCSV).split('_')[2]=='hants' or os.path.basename(outCSV).split('_')[2]=='whittaker' or os.path.basename(outCSV).split('_')[2]=='savgol') :
    #            print (Date)
                # Create in memory Raster for current date
                dest = rmem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
                dest.SetGeoTransform(rGeoT)
                dest.SetProjection(Proj.ExportToWkt())
                dest.GetRasterBand(1).WriteArray(v_values[counter,:,:])
                dest.GetRasterBand(1).SetNoDataValue(np.nan)
    
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
                masked = np.ma.MaskedArray(valuesArray,mask=np.logical_or(valuesArray == nodata,np.logical_not(rfeatArray)*tree_mask_array))
    #            print (masked)
                
                mean = masked.mean()
                csvDict.setdefault(str(Date)+'Mean',[]).append(mean)
                if np.isnan(mean) == False : 
                    std = masked.std()
                    csvDict.setdefault(str(Date)+'Std',[]).append(std)
                else :
                    csvDict.setdefault(str(Date)+'Std',[]).append(np.nan)
                
    #            masked.max()
    #            masked.sum()
    #            masked.count()
                dest = None
                rfeat_ds = None
                vfeat_ds = None
                
            counter += 1
#    print (csvDict)
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
    viIndexName = variable.split('_')[2]
    
    Y = ncds.variables['Y'][:]
    X = ncds.variables['X'][:]
    xsize = len(X)
    ysize = len(Y)
    time = [Date for Date in ncds.variables['time'][:]]
    origin_time = [Date for Date in ncds.variables['time_origin'][:]]#ncds.variables['time_origin'][:]
#    print (time)
    Proj = osr.SpatialReference()
    Proj.ImportFromEPSG(32628)
    
    originX = min(X)-0.5*3.
    originY = max(Y)-0.5*3.
    
    rGeoT = (originX,3.,0,originY,0,-3.)
    
    date_index = time.index(20170618)
    ovalues_20170618 = np.array(ncds.variables["original_PRScor_NDVI_values"][date_index,:,:])
    ovalues_20170618 = np.where(ovalues_20170618==-9999.,np.nan,ovalues_20170618)
    ovalues_20170618 = ovalues_20170618 / 10000
    
    # Extract Original and Hants Values
    variable_values = np.array(ncds.variables[variable][:])
    v_values = np.where(variable_values==-9999.,np.nan,variable_values)
    v_values = v_values / 10000
    #outName
    outFolder = os.path.join(os.path.dirname(inFile),str(viIndexName)+'_aggregate_notree')
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    # CSV Name
    outCSV = os.path.join(outFolder,'agg_%s_%s_%s.csv'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    create_csv(outCSV, time, origin_time, xsize, ysize, originX, originY, Proj, rGeoT, v_values, ovalues_20170618)
    
    # CSVt File for Attribute Join in GIS Software
    outCSVt = os.path.join(outFolder,'agg_%s_%s_%s.csvt'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    with open(outCSVt,'w') as csvt_file:
        csvt_file.write('"String",') 
        for i in range(len(time)-1) :
            csvt_file.write('"Real",')
        csvt_file.write('"Real"')
        
    # Close netCDF dataset
    ncds = None


if __name__=="__main__":
    
# =============================================================================
#     NDVI & MSAVI2
# =============================================================================
    
    inFile = "H:/Stage2018/Process/LISSAGE/TIME_SERIES.nc"
    inVectorFile = "H:/Stage2018/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor_epure.shp"
    
    variables = ["original_PRScor_GDVI_values","original_PRS_GDVI_values",
                 "original_PRScor_CIGreen_values","original_PRS_CIGreen_values"]
    
#    variable = "whittaker_PRScor_NDVI_values"
#    aggregate (inFile, inVectorFile, variable)
    
    for variable in variables :
        print (variable)
        
        aggregate (inFile, inVectorFile, variable)