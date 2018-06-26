# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:38:31 2018

@author: je

Plot Time Series Original and Smoothed Values with Max VI value and ground thruth data as Sowing and Harvest Dates

"""
import os
import matplotlib.pyplot as plt
import seaborn as sns
import geopandas as gpd
import pandas as pd
from pprint import pprint
import calendar
from matplotlib.backends.backend_pdf import PdfPages


def plot_original_smoothing (Path, inVector, originalCSV, hantsCSV, whitCSV, savgCSV, outPDFName) :
    """
    This function read the parcels vector file (currently in geojson format), Original and Smoothed Values stored in CSV Files
    to subplot for each plot original, smoothed values and ground thruth data as Sowing and Harvest Dates and Max NDVI 
    from complete VI time serie and Planet VI time serie. 
    """
    # VI Index Name
    viIndexName = os.path.basename(originalCSV).split('_')[1]
    
    # Scan Folder to extract Sensor name for Date
    df_DateSensor = pd.DataFrame({'Date':[pd.to_datetime(os.path.basename(File).split('_')[1], format='%Y%m%d') for File in os.listdir(Path) if File.endswith('RESIZE.tif')],
                                  'Sensor':[os.path.basename(File).split('_')[2] for File in os.listdir(Path) if File.endswith('RESIZE.tif')]})
#    print (df_DateSensor)

    # Read inVector File into GeoDataframe
    inVector_Ds = gpd.read_file(inVector)
#    print (inVector_Ds)
    
    # Read CSV Files and Join its attributes to GeoDataframe 
    original_complete = pd.read_csv(originalCSV)
    o_join = inVector_Ds.merge(original_complete, on='ID')
#    print (o_join)
    hants_complete = pd.read_csv(hantsCSV)
    h_join = inVector_Ds.merge(hants_complete, on='ID')
#    print (h_join)
    whit_complete = pd.read_csv(whitCSV)
    w_join = inVector_Ds.merge(whit_complete, on='ID')
#    print (w_join)
    
    savg_complete = pd.read_csv(savgCSV)
    savg_join = inVector_Ds.merge(savg_complete, on='ID')
#    print (savg_join)
    
    # Time List
    times = pd.Series(pd.date_range(start='2017-05-08', end = '2017-11-19', freq='D'))
#    print (times)

    # Creating Plot PDF
    with PdfPages(os.path.join(Path,outPDFName),'a') as pdf:
        for index in range(inVector_Ds.shape[0]) :
            # Values
            # Complete Time Serie & Planet only
            values_o = []
            values_h = []
            values_w = []
            values_savg = []
            
            for Date in times :
    #            print (Date.strftime('%Y%m%d'))
                values_o.append(o_join[Date.strftime('%Y%m%d')+'_OMean'][index])
                values_h.append(h_join[Date.strftime('%Y%m%d')+'_SMean'][index])
                values_w.append(w_join[Date.strftime('%Y%m%d')+'_SMean'][index])
                values_savg.append(savg_join[Date.strftime('%Y%m%d')+'_SMean'][index])
                
                
            series_o = pd.Series(values_o)
            series_h = pd.Series(values_h)
            series_w = pd.Series(values_w)
            series_savg = pd.Series(values_savg)

            vi_min = min(series_o.min(), series_h.min(), series_w.min(), series_savg.min())
            vi_max = max(series_o.max(), series_h.max(), series_w.max(), series_savg.max())

            # Plot
            plt.style.use('ggplot')
            # New figure
            
            plt.figure()
            
            # Categorize Complete Time Series 
            Planet_dates = []
            Planet_values = []
            S2_dates = []
            S2_values = []
            RapidEye_dates = []
            RapidEye_values = []
            for i in range(df_DateSensor.shape[0]):
                Date = df_DateSensor['Date'][i]
                Sensor = df_DateSensor['Sensor'][i]
#                print (Date, Sensor)
                if Sensor == 'Planet' :
                    Planet_dates.append(Date) 
                    Planet_values.append(series_o[times[times==Date].index[0]])
                elif Sensor == 'RapidEye' :
                    RapidEye_dates.append(Date)
                    RapidEye_values.append(series_o[times[times==Date].index[0]])
                elif Sensor == 'S2' :
                    S2_dates.append(Date)
                    S2_values.append(series_o[times[times==Date].index[0]])
            
            Planet_df = pd.DataFrame({'Date':Planet_dates,'Values':Planet_values})
            S2_df = pd.DataFrame({'Date':S2_dates,'Values':S2_values})
            RapidEye_df = pd.DataFrame({'Date':RapidEye_dates,'Values':RapidEye_values})

            plt.plot(times, series_h, 'r-', label='HANTS')
            plt.plot(times, series_w, 'm-', label='Whittaker Smoother')
            plt.plot(times, series_savg, 'c-', label='Smoothing Average')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in Planet_df['Date']], Planet_df['Values'].values.tolist(), 'bo', label='PlanetScope')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in S2_df['Date']], S2_df['Values'].values.tolist(), 'go', label='Sentinel-2')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in RapidEye_df['Date']], RapidEye_df['Values'].values.tolist(), 'yo', label='RapidEye')
            
            plt.yticks(pd.np.arange(0,0.9,0.1))
            plt.xlim(pd.to_datetime('20170501', format='%Y%m%d'),pd.to_datetime('20171201', format='%Y%m%d'))
            plt.ylim(vi_min-0.1,vi_max+0.1)
            plt.xticks(pd.date_range(start='2017-05-01', end = '2017-12-01', freq='MS'),calendar.month_name[5:13], rotation = 15)
            plt.xlabel('Time')
            plt.ylabel('Values')
            plt.suptitle ('Mean %s'%viIndexName)
            plt.title('Plot %s - %s (%s, %s, %s)'%(inVector_Ds['ID'][index],inVector_Ds['CropSyst'][index],inVector_Ds['Crop_1'][index],inVector_Ds['Crop_2'][index],inVector_Ds['Crop_3'][index]))
            
            sowing = inVector_Ds['DateSemi'][index]
            if sowing is not None :
                sowing_date = pd.to_datetime(str(sowing[6:]+sowing[3:5]+sowing[:2]), format='%Y%m%d')
                index_sowing = times[times==sowing_date].index[0]
                sowing_vi_h = series_h[index_sowing]
                sowing_vi_w = series_w[index_sowing]
                sowing_vi_savg = series_savg[index_sowing]
                plt.plot([sowing_date,sowing_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Sowing date : %s (H %s, WS %s, SAvg %s)'%(sowing_date.strftime('%d %b'),round(sowing_vi_h,3), round(sowing_vi_w,3), round(sowing_vi_savg,3)))
    
            harvest = inVector_Ds['DateReco'][index]
            if harvest is not None :
                harvest_date = pd.to_datetime(str(harvest[6:]+harvest[3:5]+harvest[:2]), format='%Y%m%d')
                index_harvest = times[times==harvest_date].index[0]
                harvest_vi_h = series_h[index_harvest]
                harvest_vi_w = series_w[index_harvest]
                harvest_vi_savg = series_savg[index_harvest]
                plt.plot([harvest_date,harvest_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Harvest date : %s (H %s, WS %s, SAvg %s)'%(harvest_date.strftime('%d %b'), round(harvest_vi_h,3), round(harvest_vi_w,3), round(harvest_vi_savg,3)))
            
            max_vi_h = series_h.max()
            max_vi_h_index = series_h[series_h==max_vi_h].index[0]
            max_date_h = times[max_vi_h_index]
            
            max_vi_w = series_w.max()
            max_vi_w_index = series_w[series_w==max_vi_w].index[0]
            max_date_w = times[max_vi_w_index]
            
            max_vi_savg = series_savg.max()
            max_vi_savg_index = series_savg[series_savg==max_vi_savg].index[0]
            max_date_savg = times[max_vi_savg_index]
            
            plt.plot([max_date_h,max_date_h],[vi_min-0.1,vi_max+0.1], 'r--', label='H Max date : %s (%s)'%(max_date_h.strftime('%d %b'), round(max_vi_h,3)))
            plt.plot([max_date_w,max_date_w],[vi_min-0.1,vi_max+0.1], 'm--', label='WS Max date : %s (%s)'%(max_date_w.strftime('%d %b'), round(max_vi_w,3)))
            plt.plot([max_date_savg,max_date_savg],[vi_min-0.1,vi_max+0.1], 'c--', label='SAvg Max date : %s (%s)'%(max_date_savg.strftime('%d %b'), round(max_vi_savg,3)))
            plt.legend(loc=4, fancybox=True, framealpha=0.5)
    
            # Save
            pdf.savefig()

if __name__=="__main__":
    
    inVector = "D:/Stage/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.geojson"
    originalCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_original_PRScor.csv"
    hantsCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_hants_PRScor.csv"
    whitCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_whittaker_PRScor_iter.csv"
    savgCSV = "D:/Stage/Output/INDICES/NDVI/Mean_NDVI_SmoothAvg_PRScor_iter.csv"
    Path = "D:/Stage/Output/INDICES/NDVI" 
    outPDFName = "Plots_PRScor_iter.pdf"
    
    plot_original_smoothing (Path, inVector, originalCSV, hantsCSV, whitCSV, savgCSV, outPDFName)