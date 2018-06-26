# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 14:12:41 2018

@author: JE

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
#import osr
import scipy as sp
import os
import re
import xml.etree.ElementTree as ET
import json
from math import cos,pi
import sys
import subprocess
import time
import concurrent.futures

def get_datatype (dataset):
    
    # Get the type of the data
    gdal_dt = dataset.GetRasterBand(1).DataType
    if gdal_dt == gdal.GDT_Byte:
        datatype = 'uint8'
    elif gdal_dt == gdal.GDT_Int16:
        datatype = 'int16'
    elif gdal_dt == gdal.GDT_UInt16:
        datatype = 'uint16'
    elif gdal_dt == gdal.GDT_Int32:
        datatype = 'int32'
    elif gdal_dt == gdal.GDT_UInt32:
        datatype = 'uint32'
    elif gdal_dt == gdal.GDT_Float32:
        datatype = 'float32'
    elif gdal_dt == gdal.GDT_Float64:
        datatype = 'float64'
    elif gdal_dt == gdal.GDT_CInt16 or gdal_dt == gdal.GDT_CInt32 or gdal_dt == gdal.GDT_CFloat32 or gdal_dt == gdal.GDT_CFloat64 :
        datatype = 'complex64'
    else:
        print ('Data type unkown')
        exit()

    return datatype

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
#        print (dicMetadata)
        
        # Open Raster File and Convert DN to TOA Reflectance
        
        g=gdal.Open(os.path.join(inPath,inRaster),gdal.GA_ReadOnly)
        if g is None:
            print ('Can\'t open Analystic MS File. Please Check it.')
            sys.exit (1)
    #    print (g)
        geoT = g.GetGeoTransform()
        Proj = g.GetProjection()
        x_size = g.RasterXSize  # Raster xsize
        y_size = g.RasterYSize  # Raster ysize
        bandCount=g.RasterCount
        nodata = g.GetRasterBand(1).GetNoDataValue()

        toaArray=sp.empty((y_size,x_size,bandCount),dtype=get_datatype(g))
        for j in range(bandCount):
            dnArray=g.GetRasterBand(j+1).ReadAsArray()
            toaArray[:,:,j]=sp.round_(dnArray*dicMetadata[j+1]*10000)
        
        # Write TOA Reflectance Image in GeoTiff
        outName=inRaster[:-4]+'_TOA.tif'
    #    print (outName)
        outFolder=os.path.join(outPath,inRaster.split('_')[0][:4]+'_'+inRaster.split('_')[0][4:6]+'_'+inRaster.split('_')[0][6:8])
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
            if nodata != None :
                ds.GetRasterBand(k+1).SetNoDataValue(nodata)
        
        g=None
        ds=None
        driver=None
        

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
            sys.exit (1)
    #    print (g)
        geoT = g.GetGeoTransform()
        Proj = g.GetProjection()
        x_size = g.RasterXSize  # Raster xsize
        y_size = g.RasterYSize  # Raster ysize
        bandCount=g.RasterCount
        nodata = g.GetRasterBand(1).GetNoDataValue()
        
        toaArray=sp.empty((y_size,x_size,bandCount),dtype=get_datatype(g))
        for j in range(bandCount):
            dnArray=g.GetRasterBand(j+1).ReadAsArray()
           
            toaArray[:,:,j]=sp.round_(dnArray*0.01*(pi*(ESUN**ESUN)/dicEAI[j+1]*cos(SOLAR_ZENITH*pi/180))*10000)
            
        # Save TOA Reflectance Image in GeoTiff
        
        outName=inRaster[:-4]+'_TOA.tif'
    #    print (outName)
        outFolder=os.path.join(outPath,os.path.basename(inPath).split('_')[0][:4]+'_'+os.path.basename(inPath).split('_')[0][4:6]+'_'+os.path.basename(inPath).split('_')[0][6:8])
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
            if nodata != None :
                ds.GetRasterBand(k+1).SetNoDataValue(nodata)
        
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
        sys.exit (1)
#    print (ds)
    layer = ds.GetLayer()
# layer Extent Case
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
                    sys.exit (1)
            #    print (g)
                if compteur == 0 :
                    geoT = gi.GetGeoTransform()
                    originX = geoT[0]
                    originY = geoT[3] 
                    pixelWidth = geoT[1]
                    pixelHeight = geoT[5]
                    nodata = gi.GetRasterBand(1).GetNoDataValue()
                    
                    Proj = gi.GetProjection()
                    
                    xOffset = int((ulx- originX) / pixelWidth)
                    yOffset = int((uly - originY) / pixelHeight)
                    
                    xsize = int((lrx - ulx)/pixelWidth)  # Raster xsize col
                    ysize = int((uly - lry)/pixelWidth)  # Raster ysize lines
                    
                    newGeoT=(ulx,pixelWidth,geoT[2],uly,geoT[4],pixelHeight)
                    band10Array=sp.empty((ysize,xsize,len(lstBand10)),dtype=get_datatype(gi))
                    
                band10Array[:,:,compteur]=gi.GetRasterBand(1).ReadAsArray(xOffset,yOffset,xsize,ysize)
                compteur+=1
                gi=None
    
    # Resample 20m Images Bands 
    
    
    compteur=0
    for bandNumber in lstBand20 :
        for file in lstFiles :
            if re.match("(.*)_{}(.*)".format(bandNumber),file):
#                print (file)
                gi=gdal.Open(os.path.join(inPath,file),gdal.GA_ReadOnly)
                if gi is None:
                    print ('Can\'t open {}. Please Check it.'.format(file))
                    sys.exit (1)
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
                    resampleArray=sp.empty((ysize,xsize,len(lstBand20)),dtype=get_datatype(gi))
                    
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
        if nodata != None :
            ds.GetRasterBand(k+1).SetNoDataValue(nodata)
    
    ds=None
    driver=None
    

def Mos_Resize (inPath,outPath, inShpFile, sensor) :
    """
    Mosaiquer et découper selon l'extent du fichier vecteur Shape, les images PLanetScope ou RapidEye 
    correspondant à une même date.
    Input : Dossier contenant toutes les images à mosaiquer
    Output : Image mosaiquée et découpée selon l'emprise du Shapefile
    
    """
    # Get Study Zone Extent 
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(inShpFile, gdal.GA_ReadOnly)
    if ds is None:
        print ('Couldn\'t open ShapeFile. Please Check it.')
        sys.exit (1)
#    print (ds)
    layer = ds.GetLayer()

    extent = layer.GetExtent()
#    print (extent)
    ulx = extent[0]
    uly = extent[3]
    lrx = extent[1]
    lry = extent[2]
    
    ds.Destroy()
    
    lstFiles = [file for file in os.listdir(inPath) if file.endswith('.tif')]
#    print (lstFiles)
    
    g=gdal.Open(os.path.join(inPath,lstFiles[0]))
    nodata = g.GetRasterBand(1).GetNoDataValue()
    
    if sensor == 0 :
        outName = 'Planet_'+lstFiles[0].split('_')[0][:4]+'_'+lstFiles[0].split('_')[0][4:6]+'_'+lstFiles[0].split('_')[0][6:8]+'_MOS_RESIZE.tif'
        Command = ["gdalwarp",'-overwrite',"-tr", "%s"%3, "%s"%3, "-co", "COMPRESS=LZW", "-srcnodata", "%s"%nodata , "-dstnodata", "%s"%nodata, "-multi", "-te","%s" %ulx,"%s" %lry,"%s" %lrx,"%s" %uly]
    elif sensor == 1:
        outName = 'RapidEye_'+lstFiles[0].split('_')[1][:4]+'_'+lstFiles[0].split('_')[1][5:7]+'_'+lstFiles[0].split('_')[1][8:10]+'_MOS_RESIZE.tif'
        Command = ["gdalwarp",'-overwrite',"-tr", "%s"%5, "%s"%5, "-co", "COMPRESS=LZW", "-srcnodata", "%s"%nodata , "-dstnodata", "%s"%nodata, "-multi", "-te","%s" %ulx,"%s" %lry,"%s" %lrx,"%s" %uly]
        
    for file in lstFiles : 
        Command += [os.path.join(inPath,file)]
        
    Command+=[os.path.join(outPath,outName)]
#    print (Command)
   
    process=subprocess.Popen(Command, stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    out,err = process.communicate()
    print (out,err)
    
    
if  __name__=='__main__':
    '''
    ################################## PLANETSCOPE  #####################################
#                              TOA Reflectance, MOSAIC AND RESIZE
    inPath="E:/Stage2018/Brutes/RAPIDEYE"
    outPath="E:/Stage2018/Output/RAPIDEYE"
    inShpFile="E:/Stage2018/Extent_Zone_Corr.shp"
    sensor=1
    
    lstFolders = [os.path.join(inPath,folder) for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder))]
    lstSubfolders = []
    for parentFolder in lstFolders :
        lstSubfolders += [os.path.join(parentFolder,subfolder) for subfolder in os.listdir(parentFolder) if os.path.isdir(os.path.join(parentFolder,subfolder))]
#    print(lstSubfolders)
    
    import concurrent.futures
    
    NUM_WORKERS = 3

    def work(folderName):
        """
        Thread Function
        """
        dnToTOA (folderName, outPath, sensor)
        return folderName
    
    start_time = time.time()
    
    # Submit the jobs to the thread pool executor.
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
    # Map the futures returned from executor.submit() to their destination windows.
    # The _example.compute function modifies no Python objects and releases the GIL. It can execute concurrently.
        futures = {executor.submit(work, folderName) for folderName in lstSubfolders}
        concurrent.futures.wait(futures)
        
    end_time = time.time()        
         
    print("Time for proccess : %ssecs" % (end_time - start_time)) 
    '''
    
     ################################## RAPIDEYE #####################################
    inPath="E:/Stage2018/Output/PLANETSCOPE"
    outPath="E:/Stage2018/Output/PLANETSCOPE"
    inShpFile="E:/Stage2018/Extent_Zone_Corr.shp"
    sensor=0
    
    lstFolders = [os.path.join(inPath,folder) for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder))]
    
    NUM_WORKERS = 3

    def work(folderName):
        """
        Thread Function
        """
        Mos_Resize (folderName, outPath, inShpFile, sensor)
        return folderName
    
    start_time = time.time()
    
    # Submit the jobs to the thread pool executor.
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
    # Map the futures returned from executor.submit() to their destination windows.
    # The _example.compute function modifies no Python objects and releases the GIL. It can execute concurrently.
        futures = {executor.submit(work, folderName) for folderName in lstFolders}
        concurrent.futures.wait(futures)
        
    end_time = time.time()
#
#    ################################### SENTINEL-2 ######################################
    '''
    inPath="E:/Stage2018/Brutes/S2"
    outPath="E:/Stage2018/Output/SENTINEL-2"
    inShpFile="E:/Stage2018/Extent_Zone_Corr.shp"

    lstFolders = [os.path.join(inPath,folder) for folder in os.listdir(inPath) if os.path.isdir(os.path.join(inPath,folder))]
        
    import concurrent.futures
    
    NUM_WORKERS = 3

    def work(folderName):
        """
        Thread Function
        """
        Resample_Resize_Stack_S2 (folderName, outPath, inShpFile)
        return folderName
    
    start_time = time.time()
    
    # Submit the jobs to the thread pool executor.
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
    # Map the futures returned from executor.submit() to their destination windows.
    # The _example.compute function modifies no Python objects and releases the GIL. It can execute concurrently.
        futures = {executor.submit(work, folderName) for folderName in lstFolders}
        concurrent.futures.wait(futures)
        
    end_time = time.time()        
         
    print("Time for proccess : %ssecs" % (end_time - start_time))
    '''

