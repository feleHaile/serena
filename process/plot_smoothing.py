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


def plot_original_smoothing (Path, inVector, originalCSV_c, smoothingCSV_c, originalCSV_p, smoothingCSV_p):
    """
    This function read the parcels vector file (currently in geojson format), Original and Smoothed Values stored in CSV Files
    to subplot for each plot original, smoothed values and ground thruth data as Sowing and Harvest Dates and Max NDVI 
    from complete VI time serie and Planet VI time serie. 
    """
    
    # Scan Folder to extract Sensor name for Date
#    Date_Sensor = {os.path.basename(File).split('_')[1]:os.path.basename(File).split('_')[2] for File in os.listdir(Path) if File.endswith('RESIZE.tif')}
#    print (Date_Sensor)
    df_DateSensor = pd.DataFrame({'Date':[pd.to_datetime(os.path.basename(File).split('_')[1], format='%Y%m%d') for File in os.listdir(Path) if File.endswith('RESIZE.tif')],
                                  'Sensor':[os.path.basename(File).split('_')[2] for File in os.listdir(Path) if File.endswith('RESIZE.tif')]})
#    print (df_DateSensor)
    
    # Rename S2 to Sentinel-2 and Planet to PlanetScope
#    df_DateSensor['Sensor'].replace('S2','Sentinel-2') 
#    df_DateSensor['Sensor'].replace('Planet','PlanetScope')
#    print (df_DateSensor)
    

    # Read inVector File into GeoDataframe
    inVector_Ds = gpd.read_file(inVector)
#    print (inVector_Ds)
    
    # Read CSV Files and Join its attributes to GeoDataframe 
    original_complete = pd.read_csv(originalCSV_c)
    o_c_join = inVector_Ds.merge(original_complete, on='ID')
#    print (o_c_join)
    smoothing_complete = pd.read_csv(smoothingCSV_c)
    s_c_join = inVector_Ds.merge(smoothing_complete, on='ID')
#    print (s_c_join)
    original_planet = pd.read_csv(originalCSV_p)
    o_p_join = inVector_Ds.merge(original_planet, on='ID')
#    print (o_p_join)
    smoothing_planet = pd.read_csv(smoothingCSV_p)
    s_p_join = inVector_Ds.merge(smoothing_planet, on='ID')
#    print (s_p_join)
    
    # Time List
    times = pd.Series(pd.date_range(start='2017-05-08', end = '2017-11-19', freq='D'))
#    print (times)

    # Subplot for each rown in GeoDataFrames
    with PdfPages(os.path.join(Path,'Plots_Hants_Planet.pdf'),'a') as pdf:
        for index in range(inVector_Ds.shape[0]) :
            # Values
            # Complete Time Serie & Planet only
            values_o_c = []
            values_s_c = []
            values_o_p = []
            values_s_p = []
            for Date in times :
    #            print (Date.strftime('%Y%m%d'))
                values_o_c.append(o_c_join[Date.strftime('%Y%m%d')+'_OMean'][index])
                values_s_c.append(s_c_join[Date.strftime('%Y%m%d')+'_SMean'][index])
                values_o_p.append(o_p_join[Date.strftime('%Y%m%d')+'_OMean'][index])
                values_s_p.append(s_p_join[Date.strftime('%Y%m%d')+'_SMean'][index])
                
            series_o_c = pd.Series(values_o_c)
            series_s_c = pd.Series(values_s_c)
            series_o_p = pd.Series(values_o_p)
            series_s_p = pd.Series(values_s_p)
            
            vi_min = min(series_o_c.min(), series_s_c.min(), series_o_p.min(), series_s_p.min())
            vi_max = max(series_o_c.max(), series_s_c.max(), series_o_p.max(), series_s_p.max())
    #        value = min(series_o_c.min(), series_s_c.min(), series_o_p.min(), series_s_p.min())
    #        lst.append(value)
    #    pprint (lst)
    #    pprint (min(lst))
            
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
    #            print (Date, Sensor)
                if Sensor == 'Planet' :
                    Planet_dates.append(Date) 
                    Planet_values.append(series_o_c[times[times==Date].index[0]])
                elif Sensor == 'RapidEye' :
                    RapidEye_dates.append(Date)
                    RapidEye_values.append(series_o_c[times[times==Date].index[0]])
                elif Sensor == 'S2' :
                    S2_dates.append(Date)
                    S2_values.append(series_o_c[times[times==Date].index[0]])
            
            Planet_df = pd.DataFrame({'Date':Planet_dates,'Values':Planet_values})
            S2_df = pd.DataFrame({'Date':S2_dates,'Values':S2_values})
            RapidEye_df = pd.DataFrame({'Date':RapidEye_dates,'Values':RapidEye_values})
            
    #        print ( Planet_df['Date'])
    #        plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in Planet_df['Date']],Planet_df['Values'].values.tolist())
            
    #        print (Planet_df,S2_df,RapidEye_df)
            
            
    #        for Date in times :
    #            if df_DateSensor[df_DateSensor['Date'== Date]] :#and df_DateSensor['Sensor']==Sentinel-2
    #                print (Date)
            # Plot
            plt.style.use('ggplot')
            # New figure
            """
            plt.figure()
            
            # Subplot Complete TS
#            plt.subplot(121)
            plt.plot(times, series_s_c, 'r-', label='HANTS')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in Planet_df['Date']], Planet_df['Values'].values.tolist(), 'bo', label='PlanetScope')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in S2_df['Date']], S2_df['Values'].values.tolist(), 'go', label='Sentinel-2')
            plt.plot([pd.to_datetime(Date, format='%Y-%m-%d') for Date in RapidEye_df['Date']], RapidEye_df['Values'].values.tolist(), 'yo', label='RapidEye')
            
            plt.yticks(pd.np.arange(0,0.9,0.1))
            plt.xlim(pd.to_datetime('20170501', format='%Y%m%d'),pd.to_datetime('20171201', format='%Y%m%d'))
            plt.ylim(vi_min-0.1,vi_max+0.1)
            plt.xticks(pd.date_range(start='2017-05-01', end = '2017-12-01', freq='MS'),calendar.month_name[5:13], rotation = 15)
            plt.xlabel('Time')
            plt.ylabel('Values')
            plt.title('Mean NDVI Complete TS - Parcel %s'%inVector_Ds['ID'][index])
            
            sowing = inVector_Ds['DateSemi'][index]
            if sowing is not None :
                sowing_date = pd.to_datetime(str(sowing[6:]+sowing[3:5]+sowing[:2]), format='%Y%m%d')
                index_sowing = times[times==sowing_date].index[0]
                sowing_vi = series_s_c[index_sowing]
            harvest = inVector_Ds['DateReco'][index]
            if harvest is not None :
                harvest_date = pd.to_datetime(str(harvest[6:]+harvest[3:5]+harvest[:2]), format='%Y%m%d')
                index_harvest = times[times==harvest_date].index[0]
                harvest_vi = series_s_c[index_harvest]
            max_vi = series_s_c.max()
            max_vi_index = series_s_c[series_s_c==max_vi].index[0]
            max_date = times[max_vi_index]
    #        textstr = 'Sowing : %s (%s) \nHarvest : %s (%s) \nMax = %s '%(sowing_date.strftime('%d %b'),round(sowing_vi,3), harvest_date.strftime('%d %b'), round(harvest_vi,3),round(max_vi,3))
            plt.plot([sowing_date,sowing_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Sowing date : %s (%s)'%(sowing_date.strftime('%d %b'),round(sowing_vi,3)))
            plt.plot([harvest_date,harvest_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Harvest date : %s (%s)'%(harvest_date.strftime('%d %b'), round(harvest_vi,3)))
            plt.plot([max_date,max_date],[vi_min-0.1,vi_max+0.1], 'c--', label='Max = %s'%round(max_vi,3))
            plt.legend(loc=4, fancybox=True, framealpha=0.5)
    #        print (textstr)
    #        plt.text(0.05, 0.95, textstr, transform=plt.subplot(121).transAxes, fontsize=12, verticalalignment='top',bbox=dict(boxstyle='round', facecolor='white', alpha=0.2))
    #        plt.annotate('Sowing : %s (%s)'%(sowing_date.strftime('%d %b'),round(sowing_vi,3)), xy=(sowing_date, sowing_vi), xytext=(sowing_date,  -sowing_vi+0.5), horizontalalignment='left', verticalalignment='top')#, arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
    #        plt.annotate('Harvest : %s (%s)'%(harvest_date.strftime('%d %b'), round(harvest_vi,3)), xy=(harvest_date, harvest_vi), xytext=(harvest_date, harvest_vi+0.05),  arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
    #        plt.annotate('Max = %s'%round(max_vi,3), xy=(times[max_vi_index], max_vi), xytext=(times[max_vi_index]-pd.Timedelta('10 days'), max_vi+0.08), arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
            # Save
            pdf.savefig()
    
    
            """
            # Subplot Planet TS
#            plt.subplot(122)
            plt.figure()
            plt.plot(times, series_s_p, 'r-', label='HANTS')
            plt.plot(times, series_o_p, 'bo', label='PlanetScope')
            plt.yticks(pd.np.arange(0,0.9,0.1))
            plt.xlim(pd.to_datetime('20170501', format='%Y%m%d'),pd.to_datetime('20171201', format='%Y%m%d'))
            plt.ylim(vi_min-0.1,vi_max+0.1)
            plt.xticks(pd.date_range(start='2017-05-01', end = '2017-12-01', freq='MS'),calendar.month_name[5:13],rotation = 15)
            plt.xlabel('Time')
            plt.ylabel('Values')
            plt.title('Mean NDVI PlanetScope TS - Parcel %s'%inVector_Ds['ID'][index])
            
            sowing = inVector_Ds['DateSemi'][index]
            if sowing is not None :
                sowing_date = pd.to_datetime(str(sowing[6:]+sowing[3:5]+sowing[:2]), format='%Y%m%d')
                index_sowing = times[times==sowing_date].index[0]
                sowing_vi = series_s_p[index_sowing]
            harvest = inVector_Ds['DateReco'][index]
            if harvest is not None :
                harvest_date = pd.to_datetime(str(harvest[6:]+harvest[3:5]+harvest[:2]), format='%Y%m%d')
                index_harvest = times[times==harvest_date].index[0]
                harvest_vi = series_s_p[index_harvest]
            max_vi = series_s_p.max()
            max_vi_index = series_s_c[series_s_p==max_vi].index[0]
            max_date = times[max_vi_index]
    #        textstr = 'Sowing : %s (%s) \nHarvest : %s (%s) \nMax = %s '%(sowing_date.strftime('%d %b'),round(sowing_vi,3), harvest_date.strftime('%d %b'), round(harvest_vi,3),round(max_vi,3))
            plt.plot([sowing_date,sowing_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Sowing date : %s (%s)'%(sowing_date.strftime('%d %b'),round(sowing_vi,3)))
            plt.plot([harvest_date,harvest_date], [vi_min-0.1,vi_max+0.1], 'k--', label='Harvest date : %s (%s)'%(harvest_date.strftime('%d %b'), round(harvest_vi,3)))
            plt.plot([max_date,max_date],[vi_min-0.1,vi_max+0.1], 'c--', label='Max = %s'%round(max_vi,3))
            plt.legend(loc=4, fancybox=True, framealpha=0.5)
    #        plt.annotate('Sowing : %s (%s)'%(sowing_date.strftime('%d %b'),round(sowing_vi,3)), xy=(sowing_date, sowing_vi), xytext=(sowing_date,  -sowing_vi+0.5), arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
    #        plt.annotate('Harvest : %s (%s)'%(harvest_date.strftime('%d %b'), round(harvest_vi,3)), xy=(harvest_date, harvest_vi), xytext=(harvest_date, harvest_vi+0.05),  arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
    #        plt.annotate('Max = %s'%round(max_vi,3), xy=(times[max_vi_index], max_vi), xytext=(times[max_vi_index]-pd.Timedelta('20 days'), max_vi+0.08), arrowprops={'arrowstyle': '-|>','ls': 'dashed','color': 'black','lw':0.8}, va='center')
    #        # Show
            plt.show()
            # Save
            pdf.savefig()
                    
    #        fig, axs = plt.subplots(1, 2, figsize=(5, 5))
    ##        print (axs.shape)
    #        axs[0].plot(times, values_s_c, 'r-', label='HANTS')
    #        axs[0].plot(times, values_o_c, 'b.', label='Original data')
    ##        axs[0].yticks = pd.np.arange(0,0.8,0.1)
    #        axs[1].plot(times, values_s_p, 'r-', label='HANTS')
    #        axs[1].plot(times, values_o_p, 'c^', label='Original data')
    #        axs[0].set_yticks(pd.np.arange(0,0.9,0.1))
    ##        axs[0].set_xticks(calendar.month_name[1:13])#, rotation=45)
    #        axs[1].set_yticks(pd.np.arange(0,0.9,0.1))
    #        axs[1].set_xticks(calendar.month_name[1:13])#, rotation=45)


    
if __name__=="__main__":
    
    inVector = "G:/NDVI/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOIN.geojson"
    originalCSV_c = "G:/NDVI/Mean_NDVI_original_complete.csv"
    originalCSV_p = "G:/NDVI/Mean_NDVI_original_planet.csv"
    smoothingCSV_c = "G:/NDVI/Mean_NDVI_hants_complete.csv"
    smoothingCSV_p = "G:/NDVI/Mean_NDVI_hants_planet.csv"
    Path = "G:/NDVI" 
    plot_original_smoothing (Path, inVector, originalCSV_c, smoothingCSV_c, originalCSV_p, smoothingCSV_p)