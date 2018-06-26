library(readr)
library(ggplot2)
library(ggthemes)
library(scales)

incsv = read.csv("F:/to_plot/Plot_CSV/agg_NDVI_original_PRS_toggplot.csv",check.names=FALSE)

attach(incsv)
df = incsv[which(ID=="MAEE3"),c("Date", "Sensor")]
df = rbind(df,c("2017-07-27","RapidEye"),c("2017-07-27","Sentinel-2"),c("2017-08-06","Sentinel-2"))

p = ggplot(data=df, aes(x=as.Date(Date),y=factor(Sensor))) + geom_point(aes(shape=Sensor),size=3, color = "deeppink1")+ 
  theme_minimal() + theme(axis.text.y = element_text(face="bold"),legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_x_date(name="Date", breaks = date_breaks("months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-12-01"))) + 
  theme(panel.grid.minor.x=element_blank())

ggsave("C:/Users/User/Desktop/serena/report/materiels_methodes/chronologie.png",width=8, height=1.5)