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
#import pandas as pd
import os
import numpy as np
from netCDF4 import Dataset
import sys
#from pprint import pprint
import csv

def create_csv(outCSV, time, time2, xsize, ysize, originX, originY, Proj, rGeoT, v_values) :
    """
    Aggregate Parcels Mean Values and Write them in CSV
    
    """
    # Create Lists for Attributes Table to write GeoJSON from GeoPandas Geodataframe
#    lstGeom = []
#    original_ID = []
#    original_Projet = []
#    original_CropSyst = []
#    original_joint_SC = []
#    original_FieldType = []
#    original_joint_Semi = []
#    original_joint_Reco = []
#    original_joint_Cult = []
#    original_joint_Biom = []
#    original_joint_Rdt = []
#    original_Crop_1 = []
#    original_Crop_2 = []
#    original_Crop_3 = []
#    original_Harvest = []
#    original_Residus = []
#    original_Nb_Trees = []
#    original_Nb_Species = []
#    original_Specie_1 = []
#    original_Specie_2 = []
#    original_Specie_3 = []
#    original_Nb_Spe_1 = []
#    original_Nb_Spe_2 = []
#    original_Nb_Spe_3 = []
#    original_Terrain = []
    
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
    
    # Create CSV File
    csvfile = open(outCSV,'w')
    writer = csv.writer(csvfile, delimiter=',')
    # Write header
    writer.writerow(time2)
    
    # For each Feature, Get Geometry, attributes and Compute Zonal Stats
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        print ("Feature %s ID %s"%(i,feature.GetField('ID')))
        geometry = feature.GetGeometryRef()
        extent = geometry.GetEnvelope()
        
        # Load Feature Attribute to List
#        lstGeom.append(str(geometry))
#        original_ID.append(feature.GetField('ID'))
#        original_Projet.append(feature.GetField('Projet'))
#        original_CropSyst.append(feature.GetField('CropSyst'))
#        original_joint_SC.append(feature.GetField('joint_SC'))
#        original_FieldType.append(feature.GetField('FieldType'))
#        original_joint_Semi.append(feature.GetField('joint_Semi'))
#        original_joint_Reco.append(feature.GetField('joint_Reco'))
#        original_joint_Cult.append(feature.GetField('joint_Cult'))
#        original_joint_Biom.append(feature.GetField('joint_Biom'))
#        original_joint_Rdt.append(feature.GetField('joint_Rdt'))
#        original_Crop_1.append(feature.GetField('Crop_1'))
#        original_Crop_2.append(feature.GetField('Crop_2'))
#        original_Crop_3.append(feature.GetField('Crop_3'))
#        original_Harvest.append(feature.GetField('Harvest'))
#        original_Residus.append(feature.GetField('Residus'))
#        original_Nb_Trees.append(feature.GetField('Nb_Trees'))
#        original_Nb_Species.append(feature.GetField('Nb_Species'))
#        original_Specie_1.append(feature.GetField('Specie_1'))
#        original_Specie_2.append(feature.GetField('Specie_2'))
#        original_Specie_3.append(feature.GetField('Specie_3'))
#        original_Nb_Spe_1.append(feature.GetField('Nb_Spe_1'))
#        original_Nb_Spe_2.append(feature.GetField('Nb_Spe_2'))
#        original_Nb_Spe_3.append(feature.GetField('Nb_Spe_3'))
#        original_Terrain.append(feature.GetField('Date'))
        
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
        
        lstMean = []
        counter = 0
        for Date in time : 
#            print (Date)
            # Create in memory Raster for current date
#            if counter == 0 :
#                drv = gdal.GetDriverByName('GTiff')
#                dest0 = drv.Create('/home/je/Bureau/test.tif', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
#                dest0.SetGeoTransform(rGeoT)
#                dest0.SetProjection(Proj.ExportToWkt())
#                dest0.GetRasterBand(1).WriteArray(v_values[counter,:,:])
                
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
            lstMean.append(mean)
#            feature_stats = {
#            'min': float(masked.min()),
#            'mean': float(masked.mean()),
#            'max': float(masked.max()),
#            'std': float(masked.std()),
#            'sum': float(masked.sum()),
#            'count': int(masked.count()),
#            'fid': int(feature.GetFID())}
#            print (feature_stats)
            
            rfeat_ds = None
            vfeat_ds = None
            
            counter += 1
        
        lstMean.insert(0,feature.GetField('ID'))
        writer.writerow(lstMean)
    
    csvfile.close()
    
#    csv_df = pd.read_csv(outCSV_original)
    
#    data = {'ID':pd.Series(original_ID),'Projet':pd.Series(original_Projet),'CropSyst':pd.Series(original_CropSyst),'SC':pd.Series(original_joint_SC),
#            'FieldType':pd.Series(original_FieldType), 'DateSemi':pd.Series(original_joint_Semi), 'DateReco':pd.Series(original_joint_Reco),
#            'Culture':pd.Series(original_joint_Cult), 'Biom':pd.Series(original_joint_Biom), 'Rdt':pd.Series(original_joint_Rdt),
#            'Crop_1':pd.Series(original_Crop_1), 'Crop_2':pd.Series(original_Crop_2), 'Crop_3':pd.Series(original_Crop_3),
#            'Harvest':pd.Series(original_Harvest),'Residus':pd.Series(original_Residus), 'Nb_Trees':pd.Series(original_Nb_Trees),
#            'Nb_Species':pd.Series(original_Nb_Species), 'Specie_1':pd.Series(original_Specie_1), 'Specie_2':pd.Series(original_Specie_2),
#            'Specie_3':pd.Series(original_Specie_3), 'Nb_Spe_1':pd.Series(original_Nb_Spe_1), 'Nb_Spe_2':pd.Series(original_Nb_Spe_2),
#            'Nb_Spe_3':pd.Series(original_Nb_Spe_3),'DateTerrain':pd.Series(original_Terrain)}
    
#    data.update (csv_df.to_dict())
#    print (data)
    
#    geom_column = [loads(geom) for geom in lstGeom]
#    gdf = gpd.GeoDataFrame(data, geometry = geom_column)
#    gdf.crs = {'init': 'epsg:32628'}
#    gdf.to_file(os.path.join(os.path.dirname(inFile),'SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.geojson'), driver='GeoJSON')

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
        time2 = [str(Date)+'_OMean' for Date in time]
        time2.insert(0,"ID")
    elif variable.startswith('hants'):
        time2 = [str(Date)+'_HMean' for Date in time]
        time2.insert(0,"ID")
    elif variable.startswith('savgol'):
        time2 = [str(Date)+'_SMean' for Date in time]
        time2.insert(0,"ID")
    elif variable.startswith('whittaker'):
        time2 = [str(Date)+'_WMean' for Date in time]
        time2.insert(0,"ID")
#    print (time2)
     
    # Extract Original and Hants Values
    variable_values = np.array(ncds.variables[variable%viIndexName][:])
    v_values = np.where(variable_values==-9999.,np.nan,variable_values)
    v_values = v_values / 10000
    # CSV Name
    outCSV = os.path.join(os.path.dirname(inFile),'Mean_%s_%s_%s.csv'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    create_csv(outCSV, time, time2, xsize, ysize, originX, originY, Proj, rGeoT, v_values)
    
    # CSVt File for Attribute Join in GIS Software
    outCSVt = os.path.join(os.path.dirname(inFile),'Mean_%s_%s_%s.csvt'%(viIndexName,variable.split('_')[0],variable.split('_')[1]))
    with open(outCSVt,'w') as csvt_file:
        csvt_file.write('"String",') 
        for i in range(len(time)-1) :
            csvt_file.write('"Real",')
        csvt_file.write('"Real"')
        
    # Close netCDF dataset
    ncds = None
    
if __name__=="__main__":
    
# =============================================================================
#     NDVI
# =============================================================================
    
    inFile = "D:/Stage/Output/INDICES/NDVI/NDVI_TIME_SERIES.nc"
    inVectorFile = "D:/Stage/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"

    variables = ["original_PRScor_%s_values","hants_PRScor_%s_values","whittaker_PRScor_%s_values","savgol_PRScor_%s_values"]
    for variable in variables :
        print (variable%os.path.basename(inFile).split('_')[0])
        aggregate (inFile, inVectorFile, variable)
        
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
