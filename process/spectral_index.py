#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 11:18:16 2018

@author: je

Calcul d'indices spectraux à partir des images Planetscope, RapidEye et Sentinel-2 prétraitées 

Ordre des bandes des imageries

PlanetScope : 1-Blue 2-Green 3-Red 4-NIR
RapidEye :    1-Blue 2-Green 3-Red 4-Red Edge 5-NIR
Sentinel-2 :  1-Blue(B2) 2-Green(B3) 3-Red(B4) 4-Red Edge(B5) 5-Red Edge(B6) 6-Red Edge(B7) 7-NIR(B8) 8-Red EdgeB8A 
              9-SWIR(B11) 10-SWIR(B12)

"""

import gdal
import os
from osgeo import gdal_array
from pprint import pprint

def compute_index (inPath, outPath, lstIndex) :
    """
    Function to compute spectral index from Satellites Images (PlanetScope, RapidEye and Sentinel-2) like  :
        NDVI (Normalized Difference Vegetation Index)           -->  NIR - R / NIR + R
        CCCI (Canopy Chlorophylle Content Index)                -->  (NIR - RE / NIR + RE) / (NIR - R / NIR + R)
        CVI (Chlorophylle Vegetation Index)                     -->  NIR*(R/B²)
        CIGreen (Chlorophylle Index Green)                      -->  NIR/R -1
        ChlRE (Chlorophylle RE)                                 -->  RE / PIR
        GDVI (Green Difference Vegetation Index)                -->
        GLI (Green Leaf Index)                                  -->
        MSAVI2 (Modified Soil Adjusted Vegetation Index)        -->
        NDRE (Normalized Difference Red Edge)                   -->
        PSRI-NIR (Plant Senescence Reflectance Index -NIR)      -->
        RE (Red Edge)                                           -->
        REIP (Red Edge Inflexion Point)                         -->
        SR (Simple Ratio)                                       -->
        GVMI (Global Vegetation Moisture Index)                 -->
        NDWI (Normalized Difference Water Index)                -->
        RDI (Ratio Drought Index)                               -->
        
        EVI (Enhanced Vegetation Index)                         -->
        SAVI (Soil Adjusted Vegetation Index)                   -->
    
    """
    # Retrieve sensor name and get band data
    sensorName = os.path.basename(inPath).split('_')[0]
#    pprint (sensorName)
    
    if sensorName == "Planet" :
        
        ds = gdal.Open(inPath, gdal.GA_ReadOnly)
        
        geoT = ds.GetGeoTransform()
        Proj = ds.GetProjection()
        xsize = ds.RasterXSize  # Raster xsize
        ysize = ds.RasterYSize  # Raster ysize
        
        b = ds.GetRasterBand(1).ReadAsArray()
        g = ds.GetRasterBand(2).ReadAsArray()
        r = ds.GetRasterBand(3).ReadAsArray()
        nir = ds.GetRasterBand(4).ReadAsArray()
        
        ds = None
        
    elif sensorName == "RapidEye" :
        
        ds = gdal.Open(inPath, gdal.GA_ReadOnly)
        
        geoT = ds.GetGeoTransform()
        Proj = ds.GetProjection()
        xsize = ds.RasterXSize  # Raster xsize
        ysize = ds.RasterYSize  # Raster ysize
        
        b = ds.GetRasterBand(1).ReadAsArray()
        g = ds.GetRasterBand(2).ReadAsArray()
        r = ds.GetRasterBand(3).ReadAsArray()
        re = ds.GetRasterBand(4).ReadAsArray()
        nir = ds.GetRasterBand(5).ReadAsArray()
        
        ds = None
    
    elif sensorName == "S2" :
        
        ds = gdal.Open(inPath, gdal.GA_ReadOnly)
        
        geoT = ds.GetGeoTransform()
        Proj = ds.GetProjection()
        xsize = ds.RasterXSize  # Raster xsize
        ysize = ds.RasterYSize  # Raster ysize
        
        b = ds.GetRasterBand(1).ReadAsArray()
        g = ds.GetRasterBand(2).ReadAsArray()
        r = ds.GetRasterBand(3).ReadAsArray()
        re1 = ds.GetRasterBand(4).ReadAsArray()
        re2 = ds.GetRasterBand(5).ReadAsArray()
        re3 = ds.GetRasterBand(6).ReadAsArray()
        nir = ds.GetRasterBand(7).ReadAsArray()
        re4 = ds.GetRasterBand(8).ReadAsArray()
        swir1 = ds.GetRasterBand(9).ReadAsArray()
        swir2 = ds.GetRasterBand(10).ReadAsArray()
 
        ds = None
    
    # Index Computation
    outName = '{}_'+os.path.basename(inPath).split('_')[1]+os.path.basename(inPath).split('_')[2]+os.path.basename(inPath).split('_')[3]+'_'+os.path.basename(inPath).split('_')[0]+'.tif'
    
    if 'NDVI' in lstIndex :
        indexName = "NDVI"
        ndvi = (nir - r) / (nir + r)
#        pprint (ndvi)
        outFolder = os.path.join(outPath,indexName)
        
        if not os.path.isdir(outFolder) :
            os.makedirs(outFolder)
            
        outFile = os.path.join(outFolder,outName.format(indexName))

        dataType=gdal_array.NumericTypeCodeToGDALTypeCode(ndvi.dtype)
        save_geotiff (ndvi, outFile, xsize, ysize, geoT, Proj, dataType)

    

def save_geotiff (indexArray, outFile, xsize, ysize, geoT, Proj, dataType) :
    """
    Function to save index array in GeoTiff Images
    
    """
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(outFile, xsize, ysize, 1, dataType)
    ds.SetGeoTransform(geoT)
    ds.SetProjection(Proj)
    ds.GetRasterBand(1).WriteArray(indexArray)
    
    ds=None
    driver=None
    
    
if  __name__=='__main__':
    
#    inPath = "/home/je/Bureau/Stage/Output/MOS_RESIZE"
    outPath = "/home/je/Bureau/Stage/Output/Indices"
    lstIndex = ['NDVI']
    
#    lstFolders = [folder for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder)) and not folder=="RapidEye"]
#    for foldername in lstFolders :
##        pprint (foldername)
#        lstFiles = [file for file in os.listdir(os.path.join(inPath,foldername)) if file.endswith('.tif')]
#        for file in lstFiles :        
#            compute_index (os.path.join(inPath,foldername,file), outPath, lstIndex)
            
    inPath= "/home/je/Bureau/Stage/Output/MOS_RESIZE/PlanetScope/Planet_2017_06_04_MOS_RESIZE.tif"
    compute_index (inPath, outPath, lstIndex)
    
    
    
    

