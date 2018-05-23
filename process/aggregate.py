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
import geopandas as gpd
from shapely.wkt import loads
import pandas as pd
#import os
import numpy as np
from netCDF4 import Dataset
import sys
#from pprint import pprint
import csv

def aggregate (inFile, inVectorFile, outCSV_original, outCSV_smoothing, task) :
    """
    Agregate inFile Raster Pixels Values in Parcel (Vector File) Mean Value
    Based on perrygeo/zonal_stats.py : https://gist.github.com/perrygeo/5667173
    and gdal/ogr cookbook : https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#calculate-zonal-statistics
    Raster and Vector Files must have Same SRS
    
    """
    # Read NetCDF File 
    # Extract some parameters for Zonal Stats
    ncds = Dataset(inFile)
    
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
    
    # Extract Original and Hants Values
    original_values = ncds.variables['original_values'][:].data
    original = np.where(original_values==-9999.,np.nan,original_values)
    
    hants_values = ncds.variables['hants_values'][:].data
    hants = np.where(hants_values==-9999.,np.nan,hants_values)
#    print (original, hants)
    
    # For Each Feature, compute mean values  (Zonal Stats) by Date both for Orignal Values and Smoothing Values
    
# =============================================================================
#    Original Values
# =============================================================================
    if 'original' in task :
    
    #    # Create Lists for Attributes Table to write GeoJSON from GeoPandas Geodataframe
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
        time2 = [str(Date)+'_OMean' for Date in time]
        time2.insert(0,"ID")
    #    print (time2)
        csvfile_o = open(outCSV_original,'w')
        writer_o = csv.writer(csvfile_o, delimiter=',')
        # Write header
        writer_o.writerow(time2)
        
        # For each Feature, Get Geometry, attributes and Compute Zonal Stats
        for i in range(layer.GetFeatureCount()):
            # Feature
            feature = layer.GetFeature(i)
    #        print ("Feature %s ID %s"%(i,feature.GetField('ID')))
    #        featureDefn = layer.GetLayerDefn()
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
            
            rGeoT = (ulx,3.,0,uly,0,-3.)
            cellsize = 3.
            
            xOffset = int((ulx - originX) / cellsize)
            yOffset = int((uly - originY) / -cellsize)
            
    #        xOffset2 = int((lrx - originX) / cellsize) +1
    #        yOffset2 = int((lry - originY) / -cellsize) +1
            
            col_nb = int((lrx - ulx)/cellsize) # Raster xsize col
            row_nb = int((uly - lry)/cellsize)  # Raster ysize lines
            
    #        col_nb = xOffset2 - xOffset
    #        row_nb = yOffset2 - yOffset
            
            lstMean = []
            counter = 0
            for Date in time : 
    #            print (Date)
                # Create in memory Raster for current date
                dest = rmem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
                dest.SetGeoTransform(rGeoT)
                dest.SetProjection(Proj.ExportToWkt())
                dest.GetRasterBand(1).WriteArray(original[counter,:,:])
    
                valuesArray = dest.ReadAsArray(xOffset,yOffset,col_nb,row_nb)
                
                # Create Raster in Memory for feature Rasterization
                rfeat_ds = rmem_drv.Create('', col_nb, row_nb, 1, gdal.GDT_Byte)
                rfeat_ds.SetGeoTransform(rGeoT)
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
            writer_o.writerow(lstMean)
        csvfile_o.close()
        
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
    
# =============================================================================
#     Smoothing Values
# =============================================================================
    
    if 'smoothing' in task :
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
        time2 = [str(Date)+'_SMean' for Date in time]
        time2.insert(0,"ID")
    #    print (time2)
        csvfile_s = open(outCSV_smoothing,'w')
        writer_s = csv.writer(csvfile_s, delimiter=',')
        # Write header
        writer_s.writerow(time2)
        
        # For each Feature, Get Geometry, attributes and Compute Zonal Stats
        for i in range(layer.GetFeatureCount()):
            feature = layer.GetFeature(i)
    #        print ("Feature %s ID %s"%(i,feature.GetField('ID')))
            geometry = feature.GetGeometryRef()
            extent = geometry.GetEnvelope()
                  
    #        print (extent)
            ulx = extent[0]
            uly = extent[3]
            lrx = extent[1]
            lry = extent[2]
            
            rGeoT = (ulx,3.,0,uly,0,-3.)
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
                dest = rmem_drv.Create('', xsize, ysize, 1, gdal.GDT_Float64) # 7 is Float_64 Gdal datatype
                dest.SetGeoTransform(rGeoT)
                dest.SetProjection(Proj.ExportToWkt())
                dest.GetRasterBand(1).WriteArray(hants[counter,:,:])
    
                valuesArray = dest.ReadAsArray(xOffset,yOffset,col_nb,row_nb)
                
                # Create Raster in Memory for feature Rasterization
                rfeat_ds = rmem_drv.Create('', col_nb, row_nb, 1, gdal.GDT_Byte)
                rfeat_ds.SetGeoTransform(rGeoT)
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
                
                rfeat_ds = None
                vfeat_ds = None
                
                counter += 1
            
            lstMean.insert(0,feature.GetField('ID'))
            writer_s.writerow(lstMean)
        
        csvfile_s.close()
            
        rmem_drv = None
        vmem_drv = None
        shpds = None
        
    # Close netCDF dataset
    ncds = None


if __name__=="__main__":
    
    inFile = "G:/NDVI/NDVI.nc"
    inVectorFile = "F:/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"
#    inVectorFile = "H:/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.shp"
    outCSV_original = "G:/NDVI/Mean_NDVI_original_complete.csv"
    outCSV_smoothing = "G:/NDVI/Mean_NDVI_hants_complete.csv"
    aggregate (inFile, inVectorFile, outCSV_original, outCSV_smoothing, task=['original','smoothing'])

