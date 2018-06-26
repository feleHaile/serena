library(readr)
library(rgdal)

inVector = readOGR("/media/je/JE3/Gbodjo_2018/Data/Terrain/SimCo_2017_CLEAN_JOIN_COR_SOPHIE_ADAMA_32628_JOINCor.geojson")

# inCSV = read.csv("/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_original_PRS.csv", check.names = FALSE) # NDVI
# inCSV = read.csv("/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_original_PRS.csv", check.names = FALSE) # MSAVI2

# inCSV = read.csv("/media/je/JE3/to_plot/NDVI_aggregate/agg_NDVI_original_PRScor.csv", check.names = FALSE) # NDVI PRSCor
inCSV = read.csv("/media/je/JE3/to_plot/MSAVI2_aggregate/agg_MSAVI2_original_PRScor.csv", check.names = FALSE) # MSAVI2 PRSCor

inPath = "/media/je/SAUVEGARDE/Stage2018/Output/INDICES/NDVI"

# Jointure GeoJSON et valeurs csv
data_join = merge(x = inVector, y = inCSV, by = "ID", all = TRUE)

# récupérer le nom du capteur pour chaque date
files <- list.files(path=inPath, pattern="*.tif", full.names=F)
dates <- c() # Modifier Dates
sensor <-c()

for (i in 1:length(files)){
  dates <- c(dates,unlist(strsplit(files[i],split="_",fixed=TRUE))[2]) 
  sensor<- c(sensor,unlist(strsplit(files[i],split="_",fixed=TRUE))[3])
}

# Pour chaque parcelle récupérer les attributs caractéristiques (ID, Culture, Système, Semi, Récolte) et 
# créer un dataframe avec 4 colonnes (Date, Avg, Std, Sensor)

total <-data.frame(stringsAsFactors = FALSE)
attach(data_join@data)
for (j in 1:length(data_join@data$ID)){
  id = data_join@data$ID[j]
  CropSys = data_join@data$SC[j]
  Culture = data_join@data$Culture[j]
  Semi = as.Date(data_join@data$DateSemi[j], format="%d/%m/%Y")
  Recolte = as.Date(data_join@data$DateReco[j], format="%d/%m/%Y")
  
  df = data.frame(Date = as.Date(dates, format="%Y%m%d"), Sensor = sensor, stringsAsFactors = FALSE) # Enlever Sensor
  df$Sensor[df$Sensor=="S2"]<-"Sentinel-2"
  df$Sensor[df$Sensor=="Planet"]<-"PlanetScope"
  df$ID <- id 
  df$CropSys <- CropSys
  df$Culture <- Culture
  df$Semi <- as.Date(Semi,format="%d%m%Y")
  df$Recolte <- as.Date(Recolte,format="%d%m%Y")
  
  avg <- c()
  stdev <- c()
  for (k in 1:length(dates)){
    avgname <- paste0(dates[k],"Mean")
    stdname <- paste0(dates[k],"Std") 
    
    avg <- c(avg,data_join@data[which(ID==id),avgname])
    stdev <- c(stdev,data_join@data[which(ID==id),stdname])
  }
  
  df$Avg <- avg
  df$Std <- stdev
  total <- rbind(total,df)
}

# To CSV
write.csv(total,"/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_original_PRScor_toggplot.csv",row.names = FALSE)
