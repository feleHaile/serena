# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:12:00 2018
@author: je

Estimate Yields using Phenology Metrics extracted in Random Forest Model
"""

import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from patsy import dmatrices
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.model_selection import cross_val_score, cross_val_predict, RandomizedSearchCV
from sklearn.feature_selection import RFE, RFECV
from sklearn.preprocessing import PolynomialFeatures
from sklearn import metrics
from math import sqrt
from pprint import pprint

def read_inputs (inMetrics,inVector) : 
    """
    Read inputs files and return merged dataframe
    """
    # Read files
    metrics_df = pd.read_csv(inMetrics)
    gdf = gpd.read_file(inVector)
    # print (metrics_df.head(),gdf.head())

    join_df = gdf.merge(metrics_df,on='ID',how='left')
    join_df[["Biom","Rdt"]] = join_df[["Biom","Rdt"]].astype(float)
    join_df["CropSyst"] = join_df.FullSyst.apply(lambda row : row.split(" ")[1])
    # print(join_df.head())
    return join_df

def best_integrated_period (inMetrics,inVector, outCSV) :
    """
    Look for integrated periods that best explain biomass 
    and yields using a simple Ordinary Least Squares Linear Regression
    """
    join_df = read_inputs(inMetrics,inVector)

    #SystCult = ["Mixed Groundnut","Pure Millet","Mixed Millet"]
    SystCult = ["Groundnut","Millet"]
    Target = ["Biom","Rdt"]
    indices = ["NDVI","GDVI","CIGreen"]
    
    corDict = {}
    for i in range (len(SystCult)) :
        # df = join_df[join_df["FullSyst"]==SystCult[i]]
        df = join_df[join_df["CropSyst"]==SystCult[i]]
        for m in range(len(Target)):
            y = np.array(df[Target[m]].astype(float))
            
            for o in range (len(indices)):
                
                for j in range(20): 
                    for k in range(j+1,j+21):
                        if 5*k <= 100 :
                            X = np.array(pd.DataFrame(data=df["CUM_%s_%s_%s"%(indices[o],5*j,5*k)],columns=["CUM_%s_%s_%s"%(indices[o],5*j,5*k)]))
    
                            lm = linear_model.LinearRegression()
                            lm.fit(X,y)
                            y_predict = lm.predict(X)
    
                            corDict.setdefault('FullSyst',[]).append(SystCult[i])
                            corDict.setdefault('id_cum',[]).append("CUM_%s_%s_%s"%(indices[o],5*j,5*k))
                            corDict.setdefault("Explained_Variance",[]).append(round(metrics.explained_variance_score(y,y_predict),4))
                            corDict.setdefault('R2',[]).append(round(metrics.r2_score(y,y_predict),4))
                            corDict.setdefault('MAE',[]).append(round(metrics.mean_absolute_error(y,y_predict),4))
                            corDict.setdefault('RMSE',[]).append(round(sqrt(metrics.mean_squared_error(y,y_predict)),4))
                            corDict.setdefault('Target',[]).append(Target[m])
    
    cordf = pd.DataFrame.from_dict(corDict)
    cordf.to_csv(outCSV,index=False)

def features_interaction(inMetrics,inVector, outCSV):
    """
    Look for interaction in phenology metrics variables combination
    """
    join_df = read_inputs(inMetrics,inVector)

    # SystCult = ["Mixed Groundnut","Pure Millet","Mixed Millet"]
    SystCult = ["Groundnut","Millet"]
    Target = ["Biom","Rdt"]

    interDict = {}

    for i in range (len(SystCult)) :
        # df = join_df[join_df["FullSyst"]==SystCult[i]]
        df = join_df[join_df["CropSyst"]==SystCult[i]]
        for m in range(len(Target)):
            y = np.array(df[Target[m]].astype(float))
            
            if (SystCult[i] ==  "Groundnut" and Target[m]== "Biom"):
                x = pd.DataFrame(data=df[["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_15_100"]], columns=
                             ["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_15_100"])
                x = x.rename(columns={"CUM_GDVI_15_100": 'CUM'})
                
            elif (SystCult[i] == "Groundnut" and Target[m]== "Rdt"):
                x = pd.DataFrame(data=df[["SOS","EOS","LOS", "MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_15_20"]], columns=
                             ["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_15_20"])
                x = x.rename(columns={"CUM_GDVI_15_20": 'CUM'})
                
            elif (SystCult[i] == "Millet" and Target[m]== "Biom"):
                x = pd.DataFrame(data=df[["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_0_60"]], columns=
                             ["SOS","EOS","LOS","MAX", "AMPL", "NDVI_MAX", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_0_60"])
                x = x.rename(columns={"CUM_GDVI_0_60": 'CUM'})
                
            elif (SystCult[i] == "Millet" and Target[m]== "Rdt"):
                x = pd.DataFrame(data=df[["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_0_65"]], columns=
                             ["SOS","EOS","LOS","MAX", "NDVI_MAX", "AMPL", "RATE_SOS_MAX", "RATE_MAX_EOS", "CUM_SOS_MAX", "CUM_MAX_EOS", "CUM_GDVI_0_65"])
                x = x.rename(columns={"CUM_GDVI_0_65": 'CUM'})
                
            x_columns = list(x)
#            print (x_columns)
                
            for feature_a in x_columns:
                for feature_b in x_columns:
                    if feature_a > feature_b :
                        interDict.setdefault('FullSyst',[]).append(SystCult[i])
                        interDict.setdefault('interaction',[]).append("%s*%s"%(feature_a,feature_b))

                        x['interaction'] = x[feature_a] * x[feature_b]
                        X = np.array(pd.DataFrame(data=x["interaction"], columns=["interaction"]))

                        lm = linear_model.LinearRegression()
                        lm.fit(X,y)
                        y_predict = lm.predict(X)

                        interDict.setdefault("Explained_Variance",[]).append(round(metrics.explained_variance_score(y,y_predict),4))
                        interDict.setdefault('R2',[]).append(round(metrics.r2_score(y,y_predict),4))
                        interDict.setdefault('MAE',[]).append(round(metrics.mean_absolute_error(y,y_predict),4))
                        interDict.setdefault('RMSE',[]).append(round(sqrt(metrics.mean_squared_error(y,y_predict)),4))
                        interDict.setdefault('Target',[]).append(Target[m])

    interdf = pd.DataFrame.from_dict(interDict)
    interdf.to_csv(outCSV,index=False)
    
def check_pvalue():
    
    join_df = read_inputs(inMetrics,inVector)

    SystCult = ["Groundnut","Millet"]
    Target = ["Biom","Rdt"]
    
    for i in range (len(SystCult)) :
        df = join_df[join_df["CropSyst"]==SystCult[i]]
        for m in range(len(Target)):
            
            y = np.array(df[Target[m]].astype(float))
            
            if (SystCult[i] ==  "Groundnut" and Target[m]== "Biom"):
                # X = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_100"],columns=["EOS:CUM_GDVI_15_100"])
                X = pd.DataFrame(data=df[["EOS","CUM_GDVI_15_100"]],columns=["EOS","CUM_GDVI_15_100"])
            
            elif (SystCult[i] == "Groundnut" and Target[m]== "Rdt"):
                # X = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_20"],columns=["EOS:CUM_GDVI_15_20"])
                X = pd.DataFrame(data=df[["EOS","CUM_GDVI_15_20"]],columns=["EOS","CUM_GDVI_15_20"])
            
            elif (SystCult[i] == "Millet" and Target[m]== "Biom"):
                X = pd.DataFrame(data=df["CUM_GDVI_0_60"], columns=["CUM_GDVI_0_60"])
                
            elif (SystCult[i] == "Millet" and Target[m]== "Rdt"):
                # X = pd.DataFrame(data=df["MAX"]*df["EOS"],columns=["MAX:EOS"])
                X = pd.DataFrame(data=df[["MAX","EOS"]],columns=["MAX","EOS"])
            
            X[Target[m]] = df[Target[m]].astype(float)
            outname = Target[m]+"_"+SystCult[i]+".csv"
            X.to_csv("/home/je/Bureau/"+outname,index=False)
    
#            X = sm.add_constant(X)
#    
#            model = sm.OLS(y, X).fit()
#            print (model.summary())#, model.pvalues)
#            
#            lm = linear_model.LinearRegression()
#            lm.fit(X,y)
#            y_predict = lm.predict(X)
#            
#            print ("MAE : ",round(metrics.mean_absolute_error(y,y_predict),4))
#            print ("RMSE : ",round(sqrt(metrics.mean_squared_error(y,y_predict)),4))
#    
def plot_models():
    
    join_df = read_inputs(inMetrics,inVector)

    SystCult = ["Groundnut","Millet"]
    Target = ["Biom","Rdt"]
    
    Title2 = ["Arachide","Mil"]
    Title1 = ["Biomasses","Rendements"]
    
    for i in range (len(SystCult)) :
        df = join_df[join_df["CropSyst"]==SystCult[i]]
    
#    y = np.array(df[Target[0]].astype(float))
#    X = pd.DataFrame(data=df[Target[1]].astype(float),columns=[Target[1]])
#     
#    lm = linear_model.LinearRegression()
#    lm.fit(X,y)
#    y_predict = lm.predict(X)
#    
#    print ('R²',metrics.r2_score(y,y_predict))
#    
#    plt.figure()            
#    plt.xlabel('Biomasse végétative (kg/ha)',labelpad=5)
#    plt.ylabel('Rendement grain (kg/ha)',labelpad=5)
#    plt.scatter(y,y_predict, edgecolors=(0, 0, 0), alpha=0.8, s= 50)
#    xmin, xmax = plt.xlim()
#    ymin, ymax = plt.ylim()
#    plt.xlim([xmin,xmax])
#    plt.ylim([ymin,ymax])
#    plt.plot([xmin,xmax], [ymin,ymax], 'k--', lw=1.5)
#    
#    plt.savefig('/home/je/Bureau/Biom_vs_Rdt_Mil.png', format='png', dpi=300)
           
        
        
        for m in range(len(Target)):
            
            y = np.array(df[Target[m]].astype(float))
            
            if (SystCult[i] ==  "Groundnut" and Target[m]== "Biom"):
                # X = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_100"],columns=["EOS:CUM_GDVI_15_100"])
                X = pd.DataFrame(data=df[["EOS","CUM_GDVI_15_100"]],columns=["EOS","CUM_GDVI_15_100"])
            
            elif (SystCult[i] == "Groundnut" and Target[m]== "Rdt"):
                # X = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_20"],columns=["EOS:CUM_GDVI_15_20"])
                X = pd.DataFrame(data=df[["EOS","CUM_GDVI_15_20"]],columns=["EOS","CUM_GDVI_15_20"])
            
            elif (SystCult[i] == "Millet" and Target[m]== "Biom"):
                X = pd.DataFrame(data=df["CUM_GDVI_0_60"], columns=["CUM_GDVI_0_60"])
                
            elif (SystCult[i] == "Millet" and Target[m]== "Rdt"):
                # X = pd.DataFrame(data=df["MAX"]*df["EOS"],columns=["MAX:EOS"])
                X = pd.DataFrame(data=df[["MAX","EOS"]],columns=["MAX","EOS"])
            
            lm = linear_model.LinearRegression()
            lm.fit(X,y)
            y_predict = lm.predict(X)
            
            plt.figure()            
            plt.xlabel('Valeurs observées (kg/ha)',labelpad=5)
            plt.ylabel('Valeurs prédites (kg/ha)',labelpad=5)
            plt.scatter(y,y_predict, edgecolors=(0, 0, 0), alpha=0.8, s=50)
            xmin, xmax = plt.xlim()
            ymin, ymax = plt.ylim()
            plt.xlim([xmin,xmax])
            plt.ylim([ymin,ymax])
            plt.plot([xmin,xmax], [ymin,ymax], 'k--', lw=1.5)            
            plt.title(Title2[i])
            
            outname = Title1[m] + "_" + Title2[i]+ '.pdf'

            plt.savefig('/home/je/Bureau/'+outname, format='pdf', dpi=300)
            
    # Biomasses
#    df = join_df[join_df["CropSyst"]==SystCult[0]]
#    y1 = np.array(df[Target[0]].astype(float))
#    X1 = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_100"],columns=["EOS:CUM_GDVI_15_100"])
#    df = join_df[join_df["CropSyst"]==SystCult[1]]
#    y2 = np.array(df[Target[0]].astype(float))
#    X2 = pd.DataFrame(data=df["CUM_GDVI_0_60"], columns=["CUM_GDVI_0_60"])
#    
#    lm = linear_model.LinearRegression()
#    lm.fit(X1,y1)
#    y_predict1 = lm.predict(X1)
#    
#    lm = linear_model.LinearRegression()
#    lm.fit(X2,y2)
#    y_predict2 = lm.predict(X2)
#    
#    plt.figure()  
#    plt.subplot(121)          
#    plt.xlabel('Valeurs observées (kg/ha)',labelpad=5)
#    plt.ylabel('Valeurs prédites (kg/ha)',labelpad=5)
#    plt.scatter(y1,y_predict1, edgecolors=(0, 0, 0), alpha=0.8)
#    xmin, xmax = plt.xlim()
#    ymin, ymax = plt.ylim()
#    plt.xlim([xmin,xmax])
#    plt.ylim([ymin,ymax])
#    plt.plot([xmin,xmax], [ymin,ymax], 'k--', lw=1.5)            
#    plt.title(Title1[0] + " " + Title2[0],loc="left")
#    
#    plt.subplot(122)
#    plt.xlabel('Valeurs observées (kg/ha)',labelpad=5)
#    plt.ylabel('Valeurs prédites (kg/ha)',labelpad=5)
#    plt.scatter(y2,y_predict2, edgecolors=(0, 0, 0), alpha=0.8)
#    xmin, xmax = plt.xlim()
#    ymin, ymax = plt.ylim()
#    plt.xlim([xmin,xmax])
#    plt.ylim([ymin,ymax])
#    plt.plot([xmin,xmax], [ymin,ymax], 'k--', lw=1.5)            
#    plt.title(Title1[m] + " " + Title2[i],loc="left")
#    
#    outname = Title1[0] + " " + Title2[1]+ '.png'
#
#    plt.savefig('/home/je/Bureau/Mod_Biom.png', format='png', dpi=300)
#    
#    # Rendements 
#    y = np.array(df[Target[1]].astype(float))
#    X1 = pd.DataFrame(data=df["EOS"]*df["CUM_GDVI_15_20"],columns=["EOS:CUM_GDVI_15_20"])
#    X2 = pd.DataFrame(data=df["MAX"]*df["EOS"],columns=["MAX:EOS"])
            
            
            




if  __name__=='__main__':

    inMetrics = "/media/je/SAUVEGARDE/Cours_SIGMA/1001_Stage/Stage2018/Process/METRICS/METRICS_NDVI_GDVI_CIGreen_whittaker_PRScor.csv"
    inVector = "/media/je/SAUVEGARDE/Cours_SIGMA/1001_Stage/Stage2018/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor_epure.geojson"


    # best_integrated_period (inMetrics,inVector,outCSV="/media/je/SAUVEGARDE/Stage2018/Process/METRICS/Cumuls_NDVI_GDVI_CIGreen_UniqueMil.csv")
    # check_pvalue()
    # features_interaction(inMetrics,inVector, outCSV="/media/je/SAUVEGARDE/Stage2018/Process/METRICS/Interactions_NDVI_GDVI_CIGreen_UniqueMil.csv")
    plot_models()
    