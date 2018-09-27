library(ggplot2)
library(readr)
library(rgdal)
library(foreach)

inVector = readOGR("/media/je/JE3/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor.geojson")

fileName <- c("/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_whittaker_PRS.csv",
              "/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_whittaker_PRScor.csv")
              # "/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_hants_PRS.csv",
              # "/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_hants_PRScor.csv",
              # "/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_hants_PRS.csv",
              # "/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_hants_PRScor.csv",
              # "/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_whittaker_PRS.csv",
              # "/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_whittaker_PRScor.csv")

outFileName <- c("/media/je/JE3/to_plot/Plot_CSV/agg_NDVI_whittaker_PRS_toggplot.csv",
                 "/media/je/JE3/to_plot/Plot_CSV/agg_NDVI_whittaker_PRScor_toggplot.csv")
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_NDVI_hants_PRS_toggplot.csv",
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_NDVI_hants_PRScor_toggplot.csv",
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_hants_PRS_toggplot.csv",
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_hants_PRScor_toggplot.csv",
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_whittaker_PRS_toggplot.csv",
                 # "/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_whittaker_PRScor_toggplot.csv")

# lister les dates
dates <- seq(as.Date("2017/05/08"), by = "day", length.out = 196)

foreach(i=1:length(fileName)) %do% {
  inCSV <- read.csv(fileName[i], check.names = FALSE)
  # Jointure GeoJSON et valeurs csv
  data_join <- merge(x = inVector, y = inCSV, by = "ID", all = TRUE)
  
  # créer un dataframe avec 4 colonnes (Date, Avg, Std)
  
  total <- data.frame(stringsAsFactors = FALSE)

  for (j in 1:length(data_join@data$ID)){
    id = data_join@data$ID[j]
    
    df = data.frame(Date = as.Date(dates, format="%Y%m%d"), stringsAsFactors = FALSE) 
    
    avg <- c()
    stdev <- c()
    for (k in 1:length(dates)){
      cdate = as.character(dates[k])
      cdp=paste(substr(cdate,1,4),substr(cdate,6,7),substr(cdate,9,10), sep = "")
      avgname <- paste0(cdp,"Mean") 
      stdname <- paste0(cdp,"Std")
      
      avg <- c(avg,data_join@data[which(data_join$ID==id),avgname])
      stdev <- c(stdev,data_join@data[which(data_join$ID==id),stdname])
    }
    
    df$Avg <- avg
    df$Std <- stdev
    df$ID <- id
    total <- rbind(total,df)
  }
  
  # To CSV
  write.csv(total,outFileName[i],row.names = FALSE)
}
