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
import json
#from lxml import etree
from multiprocessing import Pool
from math import cos,pi

def dnToTOA (inPath, outPath, sensor) :
    """
    Convertir les DNs en Radiance puis réflectance des images PlanetScope et RapidEye
    sensor = 0 if PlanetScope or 1 if RapidEye
    
    """
    if sensor == 0 : # PlanetScope
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
            elif re.match("(.*)reflectanceCoefficient",element.tag):
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
            dnArray=g.GetRasterBand(j+1).ReadAsArray()
            toaArray[:,:,j]=dnArray*dicMetadata[j+1]
            toaArray[:,:,j]=sp.where(toaArray[:,:,j]==0,sp.nan,toaArray[:,:,j])
            toaArray=sp.round_(toaArray*10000)
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

    elif sensor == 1 : # RapidEye
        # Search MS Image and Metadata XML File
        lstFile= [file for file in os.listdir(inPath)]
    #    print (lstFile)
        for file in lstFile:
            if file.split('_')[-1]=="Analytic.tif":
                inRaster=file
            elif file.split('_')[-1]=="metadata.json":
                inMetadata=file
#        print (inRaster,inMetadata)
                
        # Open Raster File and Convert DN to TOA Reflectance
        
          # Exo-Atmospheric Irradiance (EAI) in W/m²/µm for each band successively 1-B 2-G 3-R 4-RE 5-NIR
        dicEAI={1:1997.8,2:1863.5,3:1560.4,4:1395.,5:1124.4}
        
          # Solar Zenithal Angle (90°-sun elevation) Sun elevation is found in Metadata json file
        with open(os.path.join(inPath,inMetadata)) as json_file :
            content=json.load(json_file)
            SOLAR_ZENITH = 90-float(content['properties']['sun_elevation'])
            ACQUISITION_DATE=content['properties']['acquired'][:10]
#            print (SOLAR_ZENITH)
#        print  (ACQUISITION_DATE)
        
          # Compute Earth-Sun Distance at the day of acquisition in Astronomical Units.
           # At this time, the Earth-Sun Distance is done on this website : http://www.instesre.org/Aerosols/sol_calc.htm
        dicESUN={"2017-07-27":1.015107862579904,"2017-10-06":0.99966670968068246,
                 "2017-10-26":0.994036745275356,"2017-11-09":0.990452817341342}
        
        ESUN=float(dicESUN[ACQUISITION_DATE])
#        print (ESUN)
        
         # Create Gdal Dataset to compute TOA Reflectance
        
        g=gdal.Open(os.path.join(inPath,inRaster),gdal.GA_ReadOnly)
        if g is None:
            print ('Can\'t open Analystic MS File. Please Check it.')
    #    print (g)
        geoT = g.GetGeoTransform()
        Proj = g.GetProjection()
        x_size = g.RasterXSize  # Raster xsize
        y_size = g.RasterYSize  # Raster ysize
        bandCount=g.RasterCount
        
        radArray=sp.empty((y_size,x_size,bandCount),dtype=sp.float64)
        toaArray=sp.empty((y_size,x_size,bandCount),dtype=sp.float64)
        for j in range(bandCount):
            dnArray=g.GetRasterBand(j+1).ReadAsArray()
            radArray[:,:,j]=dnArray*0.01
            toaArray[:,:,j]=radArray[:,:,j]*(pi*ESUN**ESUN/dicEAI[j+1]*cos(SOLAR_ZENITH*pi/180))
        toaArray=sp.where(toaArray==0,sp.nan,toaArray)
        toaArray=sp.round_(toaArray*10000)
        
        # Save TOA Reflectance Image in GeoTiff
        
        outName=inRaster[:-4]+'_TOA.tif'
    #    print (outName)
        outFolder=os.path.join(outPath,os.path.basename(inPath).split('_')[0])
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
        
        
if  __name__=='__main__':
    
    inPath="/home/je/Bureau/Stage/Gbodjo_2018/Data/Brutes/RAPIDEYE/2017_07_27"
#    outPath="/home/je/Bureau/Stage/Output/RapidEye"
#    sensor=1
#    dnToTOA(inPath, outPath, sensor)
    def process_folder(foldername) :
        outPath="/home/je/Bureau/Stage/Output/RapidEye"
        sensor=1
        
        dnToTOA(os.path.join(inPath,foldername), outPath, sensor)
        return foldername

    pool = Pool(processes=4)
    lstFolders = [folder for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder))]
#    print (lstFolders)
    results = pool.imap(process_folder, lstFolders)
   
#    for result in results :
#        print (result)