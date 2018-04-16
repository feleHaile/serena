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
import ogr
import scipy as sp
import os
import re
import xml.etree.ElementTree as ET
import json
#from lxml import etree
#from multiprocessing import Pool
from math import cos,pi

def dnToTOA (inPath, outPath, sensor) :
    """
    Convertir les DNs en Radiance et/ou réflectance des images PlanetScope et RapidEye
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
        
        g=None
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
            print ('Can\'t open Analystic File. Please Check it.')
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
        
        g=None
        ds=None
        driver=None
        

def Resample_Resize_Stack_S2 (inPath, outPath, inShpFile, option="SRE"):
    """
    Rééchantillonner les bandes Sentinel-2 de 20m à 10m de résolution spatiale par la méthode du Plus Proche Voisin 
    Découper selon l'extent du Shapefile
    et les stacker avec les bandes originelles à 10m
    Input : Dossier image contenant les différentes bandes
    Seules les bandes SRE seront utilisées
    10m -- B2, B3, B4, B8
    20m -- B5, B6, B7, B8A, B11 et B12
    
    Output : 1 Image avec toutes les bandes à  10m 
    avec cet ordre : B2, B3, B4, B5, B6, B7, B8, B8A, B11 et B12 
    
    """
    
    # Get Study Zone Extent
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(inShpFile, gdal.GA_ReadOnly)
    if ds is None:
        print ('Couldn\'t open ShapeFile. Please Check it.')
#    print (ds)
    layer = ds.GetLayer()
    extent = layer.GetExtent()
#    print (extent)
    ulx = extent[0]
    uly = extent[3]
    lrx = extent[1]
    lry = extent[2]
    
    ds.Destroy()
    
    # Define some lists
    if option == "SRE":
        lstFiles= [file for file in os.listdir(inPath) if file.endswith('.tif') and file.split('_')[-2]=="SRE"]
    elif option == "FRE":
        lstFiles= [file for file in os.listdir(inPath) if file.endswith('.tif') and file.split('_')[-2]=="FRE"]
        
#    print (lstFiles)
    lstBand10 = ["B2", "B3", "B4", "B8"]
    lstBand20 = ["B5", "B6", "B7", "B8A", "B11", "B12"]
    
    # Create an array for 10 m band data
    
    compteur=0
    for bandNumber in lstBand10 :
        for file in lstFiles :
            if re.match("(.*)_{}(.*)".format(bandNumber),file) and not file.endswith('B8A.tif'):
                gi=gdal.Open(os.path.join(inPath,file),gdal.GA_ReadOnly)
                if gi is None:
                    print ('Can\'t open {}. Please Check it.'.format(file))
            #    print (g)
                if compteur == 0 :
                    geoT = gi.GetGeoTransform()
                    originX = geoT[0]
                    originY = geoT[3] 
                    pixelWidth = geoT[1]
                    pixelHeight = geoT[5]
                    
                    Proj = gi.GetProjection()
                    
                    xOffset = int((ulx- originX) / pixelWidth)
                    yOffset = int((uly - originY) / pixelHeight)
                    
                    xsize = int((lrx - ulx)/pixelWidth)  # Raster xsize col
                    ysize = int((uly - lry)/pixelWidth)  # Raster ysize lines
                    
                    newGeoT=(ulx,pixelWidth,geoT[2],uly,geoT[4],pixelHeight)
                    band10Array=sp.empty((ysize,xsize,len(lstBand10)),dtype=sp.int16)
                    
                band10Array[:,:,compteur]=gi.GetRasterBand(1).ReadAsArray(xOffset,yOffset,xsize,ysize)
                compteur+=1
                gi=None
    
    # Resample 20m Images Bands 
    
    resampleArray=sp.empty((ysize,xsize,len(lstBand20)),dtype=sp.int16)
    compteur=0
    for bandNumber in lstBand20 :
        for file in lstFiles :
            if re.match("(.*)_{}(.*)".format(bandNumber),file):
#                print (file)
                gi=gdal.Open(os.path.join(inPath,file),gdal.GA_ReadOnly)
                if gi is None:
                    print ('Can\'t open {}. Please Check it.'.format(file))
            #    print (g)
                if compteur == 0 :
                    geoT2 = gi.GetGeoTransform()
                    originX2 = geoT2[0]
                    originY2 = geoT2[3] 
                    pixelWidth2 = geoT2[1]
                    pixelHeight2 = geoT2[5]
                    
                    xOffset2 = int((ulx- originX2) / pixelWidth2)
                    yOffset2 = int((uly - originY2) / pixelHeight2)
                    
                    xsize2 = int((lrx - ulx)/pixelWidth2)  # Raster xsize col
                    ysize2 = int((uly - lry)/pixelWidth2)  # Raster ysize lines
                    
                biArray=gi.GetRasterBand(1).ReadAsArray(xOffset2,yOffset2,xsize2,ysize2)
                dataType=gdal_array.NumericTypeCodeToGDALTypeCode(biArray.dtype)
            #    print (dataType)
                
                # Create new raster to resample		
                mem_drv = gdal.GetDriverByName('MEM')
                dest = mem_drv.Create('', xsize, ysize, 1, dataType)
                dest.SetGeoTransform(newGeoT)
                dest.SetProjection(Proj)
                
                gdal.ReprojectImage(gi, dest, Proj, Proj ,gdal.GRA_NearestNeighbour)
                
                resampleArray[:,:,compteur]=dest.GetRasterBand(1).ReadAsArray()
                compteur+=1
                gi=None
                dest=None
    
#    print (resampleArray)
    
    # Stack all Bands and save to GeoTiff images in this order : B2, B3, B4, B5, B6, B7, B8, B8A, B11 et B12
    StackArray=sp.empty((ysize,xsize,10),dtype=sp.int16)
    StackArray[:,:,:3] = band10Array[:,:,:3] # B2, B3, B4 
    StackArray[:,:,3:6] = resampleArray[:,:,:3] # B5, B6, B7
    StackArray[:,:,6] = band10Array[:,:,3] # B8
    StackArray[:,:,7:10] = resampleArray[:,:,3:6] # B8A, B11 et B12
    
    # Save Stack Image
    outName='S2'+'_'+os.path.basename(inPath).split('_')[1][:4]+'_'+os.path.basename(inPath).split('_')[1][4:6]+'_'
    outName+=os.path.basename(inPath).split('_')[1][6:8]+'_RESAMPLE_RESIZE_STACK.tif'
    
#    print (outName)
    outFolder=os.path.join(outPath)#,os.path.basename(inPath).split('_')[0])
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
  
    outFile=os.path.join(outFolder,outName)
    
    driver=gdal.GetDriverByName('GTiff')
    dataType=gdal_array.NumericTypeCodeToGDALTypeCode(StackArray.dtype)
#    dataType=gdal.GetDataTypeName(toaArray.DataType)
    ds=driver.Create(outFile,xsize,ysize,StackArray.shape[2],dataType,['COMPRESS=LZW'])
    ds.SetGeoTransform(newGeoT)
    ds.SetProjection(Proj)
    for k in range(StackArray.shape[2]):
        ds.GetRasterBand(k+1).WriteArray(StackArray[:,:,k])
    
    ds=None
    driver=None
        
if  __name__=='__main__':
    """
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
    """
    ##################################################################################
    inPath="/home/je/Bureau/Stage/Gbodjo_2018/Data/Brutes/S2/SENTINEL2A_20170727-114348-215_L2A_T28PCB_D_V1-4"
    outPath="/home/je/Bureau/Stage/Output/Sentinel-2"
    inShpFile="/home/je/Bureau/Stage/Gbodjo_2018/Data/Process/Extent_Zone.shp"
    Resample_Resize_Stack_S2 (inPath, outPath, inShpFile)