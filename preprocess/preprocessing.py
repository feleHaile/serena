#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:52:19 2018

@author: je

Fonctions de conversion de DN en Radiance puis Reflectance TOA puis Mosaiquage et Clip 
pour les images :
    -PlanetScope
    -RapidEye

Fonctions de Rééchantillonnage et Clip 
pour les images :
    -Sentinel-2

Fonctions de calcul des indices spectraux /et texturaux?

"""

import gdal
from osgeo import gdal_array
import scipy as sp
import os
import re
import xml.etree.ElementTree as ET
#from lxml import etree
from multiprocessing import Pool

def dnToTOA (inPath, outPath, sensor) :
    """
    Convertir les DNs en Radiance puis réflectance des images PlanetScope et RapidEye
    sensor = 0 if PlanetScope or 1 if RapidEye
    
    """
    # Search MS Image and Metadata XML File
    lstFile= [file for file in os.listdir(inPath)]
#    print (lstFile)
    for file in lstFile:
        if file.split('_')[-1]=="AnalyticMS.tif":
            inRaster=file
        elif file.split('_')[-1]=="metadata.xml":
            inMetadata=file
#    print (inRaster,inMetadata)
    
    # Parse XML and Get reflectanceCoefficient from Metadata File
    tree = ET.parse(os.path.join(inPath,inMetadata))
    lstBand=[]
    lstRCoef=[]
    dicMetadata={}
    for element in tree.iter():
        if re.match("(.*)bandNumber",element.tag):
            lstBand.append(int(element.text))
        if re.match("(.*)reflectanceCoefficient",element.tag):
            lstRCoef.append(float(element.text))
#    print (lstBand, lstRCoef)
    for i in range(len(lstBand)):
        dicMetadata.update({lstBand[i]:lstRCoef[i]})
#    print (dicMetadata)
    
    # Open Raster File and Convert DN to TOA Reflectance
    
    g=gdal.Open(os.path.join(inPath,inRaster),gdal.GA_ReadOnly)
    if g is None:
        print ('Can\'t open Analystic MS File. Please Check it.')
#    print (g)
    geoT = g.GetGeoTransform()
    Proj = g.GetProjection()
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize
    bandCount=g.RasterCount
    
    toaArray=sp.empty((y_size,x_size,bandCount),dtype=sp.float64)
    for j in range(bandCount):
        toaArray[:,:,j]=g.GetRasterBand(j+1).ReadAsArray()
        toaArray[:,:,j]=toaArray[:,:,j]*dicMetadata[j+1]*10000
        toaArray[:,:,j]=sp.where(toaArray[:,:,j]==0,sp.nan,toaArray[:,:,j])
        toaArray=sp.round_(toaArray)
    #toaArray=sp.where(toaArray==0,sp.nan,toaArray)
    
    # Write TOA Reflectance Image in GeoTiff
    outName=inRaster[:-4]+'_TOA.tif'
#    print (outName)
    outFolder=os.path.join(outPath,inRaster.split('_')[0])
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
        
    outFile=os.path.join(outFolder,outName)
    
    driver=gdal.GetDriverByName('GTiff')
    dataType=gdal_array.NumericTypeCodeToGDALTypeCode(toaArray.dtype)
#    dataType=gdal.GetDataTypeName(toaArray.DataType)
    ds=driver.Create(outFile,x_size,y_size,bandCount,dataType,['COMPRESS=LZW'])
    ds.SetGeoTransform(geoT)
    ds.SetProjection(Proj)
    for k in range(bandCount):
        ds.GetRasterBand(k+1).WriteArray(toaArray[:,:,k])
    ds=None
    driver=None
    
#    print (datatype)

if  __name__=='__main__':
    
    inPath="/home/je/Bureau/Stage/Gbodjo_2018/Data/Brutes/PLANET/2017_07_27"
    
    def process_folder(foldername) :
        outPath="/home/je/Bureau/Stage/Output/PlanetScope"
        sensor=0
        
        dnToTOA(os.path.join(inPath,foldername), outPath, sensor)
        return foldername

    pool = Pool(processes=4)
    lstFolders = [folder for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder))]
    
    results = pool.imap(process_folder, lstFolders)
   
    for result in results :
        print (result)