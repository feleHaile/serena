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

import rasterio
import subprocess
import os
import numpy as np
from pprint import pprint
from multiprocessing import Pool
import time

def compute_indices (inPath, outPath, lstIndex) :
    """
    Function to compute spectral index from Satellites Images (PlanetScope, RapidEye and Sentinel-2) like  :
        - NDVI (Normalized Difference Vegetation Index)           -->  nir - red / nir + red
        - CCCI (Canopy Chlorophylle Content Index)                -->  (nir - red-edge / nir + red-edge) / (nir - red / nir + red)
        - CVI (Chlorophylle Vegetation Index)                     -->  nir * (red / blue²)
        - CIGreen (Chlorophylle Index Green)                      -->  nir / red - 1
        - ChlRE (Chlorophylle RE)                                 -->  red-edge / nir
        - GDVI (Green Difference Vegetation Index)                -->
        - GLI (Green Leaf Index)                                  -->
        - MSAVI2 (Modified Soil Adjusted Vegetation Index)        -->  (2 * (nir + 1) - sqrt((2 * nir + 1)² - 8 * (nir - red)))/2
        - NDRE (Normalized Difference Red Edge)                   -->
        - PSRI-NIR (Plant Senescence Reflectance Index -NIR)      -->
        - RE (Red Edge)                                           -->
        - REIP (Red Edge Inflexion Point)                         -->
        - SR (Simple Ratio)                                       -->
        - GVMI (Global Vegetation Moisture Index)                 -->
        - NDWI (Normalized Difference Water Index)                -->
        - RDI (Ratio Drought Index)                               -->
        
        - EVI (Enhanced Vegetation Index)                         -->  
        - SAVI (Soil Adjusted Vegetation Index)                   -->
    
    """
    # Retrieve sensor name and get band data
    sensorName = os.path.basename(inPath).split('_')[0]
#    pprint (sensorName)
    
    with rasterio.open (inPath) as ds :
        
        nodataValue = ds.nodata

        if sensorName == "Planet" :
            b = ds.read(1)/10000
            g = ds.read(2)/10000
            r = ds.read(3)/10000
            nir = ds.read(4)/10000
            
            if nodataValue != None :
                
                b=np.where(b==nodataValue,np.nan,b)
                g=np.where(g==nodataValue,np.nan,g)
                r=np.where(r==nodataValue,np.nan,r)
                nir=np.where(nir==nodataValue,np.nan,nir)
                
        elif sensorName == "RapidEye" :
            
            b = ds.read(1)/10000
            g = ds.read(2)/10000
            r = ds.read(3)/10000
            re = ds.read(4)/10000
            nir = ds.read(5)/10000
            
            if nodataValue != None :
                
                b=np.where(b==nodataValue,np.nan,b)
                g=np.where(g==nodataValue,np.nan,g)
                r=np.where(r==nodataValue,np.nan,r)
                re=np.where(re==nodataValue,np.nan,re)
                nir=np.where(nir==nodataValue,np.nan,nir)

        elif sensorName == "S2" :
            
            b = ds.read(1)/10000
            g = ds.read(2)/10000
            r = ds.read(3)/10000
            re1 = ds.read(4)/10000
            re2 = ds.read(5)/10000
            re3 = ds.read(6)/10000
            nir = ds.read(7)/10000
            re4 = ds.read(8)/10000
            swir1 = ds.read(9)/10000
            swir2 = ds.read(10)/10000
            
            if nodataValue != None :
                
                b=np.where(b==nodataValue,np.nan,b)
                g=np.where(g==nodataValue,np.nan,g)
                r=np.where(r==nodataValue,np.nan,r)
                re1=np.where(re1==nodataValue,np.nan,re1)
                re2=np.where(re2==nodataValue,np.nan,re2)
                re3=np.where(re3==nodataValue,np.nan,re3)
                nir=np.where(nir==nodataValue,np.nan,nir)
                re4=np.where(re4==nodataValue,np.nan,re4)
                swir1=np.where(swir1==nodataValue,np.nan,swir1)
                swir2=np.where(swir2==nodataValue,np.nan,swir2)
                
                
        profile = ds.profile
        profile.update (dtype=rasterio.float64, count=1, nodata=np.nan, compress='lzw')
#        pprint (profile)
    
    ############ Index Computation #############
    outName1 = '{}_'+os.path.basename(inPath).split('_')[1]+os.path.basename(inPath).split('_')[2]+os.path.basename(inPath).split('_')[3]+'_'+os.path.basename(inPath).split('_')[0]+'.tif'
    outName2 = '{}_'+os.path.basename(inPath).split('_')[1]+os.path.basename(inPath).split('_')[2]+os.path.basename(inPath).split('_')[3]+'_'+os.path.basename(inPath).split('_')[0]+'_RESAMPLE_RESIZE.tif'
    
# =============================================================================
#     NDVI
# =============================================================================
    if 'NDVI' in lstIndex :
        indexName = "NDVI"
        ndvi = (nir - r) / (nir + r)
#        pprint (ndvi)
        outFolder = os.path.join(outPath,indexName)
        
        if not os.path.isdir(outFolder) :
            os.makedirs(outFolder)
            
        outFile1 = os.path.join(outFolder,outName1.format(indexName))
        outFile2 = os.path.join(outFolder,outName2.format(indexName))
        
        with rasterio.open(outFile1,'w', **profile) as outDS :
           outDS.write(ndvi.astype(rasterio.float64),1)
        
        Resample_SaveGeoTiff(outFile1,outFile2,nodata=np.nan)
    
# =============================================================================
#     MSAVI
# =============================================================================
            
    if 'MSAVI2' in lstIndex :
        indexName = "MSAVI2"
        msavi2 = (2 * (nir + 1) - np.sqrt((2 * nir + 1)**2 - 8 * (nir - r)))/2
#        pprint (msavi2)
        outFolder = os.path.join(outPath,indexName)
        if not os.path.isdir(outFolder) :
            os.makedirs(outFolder)
            
        outFile1 = os.path.join(outFolder,outName1.format(indexName))
        outFile2 = os.path.join(outFolder,outName2.format(indexName))
        
        with rasterio.open(outFile1,'w', **profile) as outDS :
           outDS.write(msavi2.astype(rasterio.float64),1)
        
        Resample_SaveGeoTiff(outFile1,outFile2,nodata=np.nan)

# =============================================================================
#      EVI
# =============================================================================
        
def Resample_SaveGeoTiff (outFile1,outFile2,nodata) :
    """
    Function to Resample and Resize indices images for all sensors to 3 meters spatial resolution.
    The upper left and lower right coordinates are derived from equations below :
        xsize = ( lrx - ulx ) / pixelWidth
        ysize = ( lry - uly) / pixelHeight
    Origin of images are fixed at (ulx,uly) and Numbers of Colums and Lines are respectively 6500 and 6460
    GdalWarp Command Line Application are used for the process.
    """
    ulx = 332673.
    uly = 1618706.
    lrx = 352173.
    lry = 1599326.
      
    Command = ["gdalwarp","-r","near", '-overwrite',"-tr", "%s"%3, "%s"%3, "-co", "COMPRESS=LZW", "-srcnodata", "%s"%nodata , "-dstnodata", "%s"%nodata, "-te","%s" %ulx,"%s" %lry,"%s" %lrx,"%s" %uly]
    Command += ["%s"%outFile1, "%s"%outFile2]
    
    process=subprocess.Popen(Command, stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    out,err = process.communicate()
    print (out,err)
    
    os.remove(outFile1)
    
if  __name__=='__main__':
    
    inPath= "/home/je/Bureau/Stage/Output/MOS_RESIZE/SENTINEL-2"
    outPath = "/home/je/Bureau/Stage/Output/INDICES"
    lstIndex = ['NDVI','MSAVI2']
    
    lstFiles = [os.path.join(inPath,file) for file in os.listdir(inPath) if file.endswith(".tif")]
#    lstFiles = ["/home/je/Bureau/Stage/Output/MOS_RESIZE/SENTINEL-2/S2_2017_07_27_RESAMPLE_RESIZE_STACK.tif",
#                "/home/je/Bureau/Stage/Output/MOS_RESIZE/RAPIDEYE/RapidEye_2017_07_27_MOS_RESIZE.tif",
#                "/home/je/Bureau/Stage/Output/MOS_RESIZE/PLANETSCOPE/Planet_2017_07_27_MOS_RESIZE.tif"]
#    pprint (lstFiles)
#    compute_indices (lstFiles[0], outPath, lstIndex)
            
    def process_file(fileName) :
        compute_indices (fileName, outPath, lstIndex)
        return fileName
    
    pool = Pool()
    
    start_time = time.time()
    
    results = pool.imap(process_file, lstFiles)
    
    for result in results :
        print (result)
    
    end_time = time.time() 
    
    pprint("Time for proccess : %ssecs" % (end_time - start_time))