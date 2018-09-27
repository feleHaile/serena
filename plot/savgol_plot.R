library(ggplot2)
library(readr)
library(foreach)
library(scales)

prs = read.csv("/media/je/JE3/R_plot/NDVI_aggregate/pxl_NDVI_savgol_PRS.csv",colClasses = 
                 c("character","double","double","numeric","numeric","integer","double","integer","integer"))

prs_cor = read.csv("/media/je/JE3/R_plot/NDVI_aggregate/pxl_NDVI_savgol_PRScor.csv",colClasses = 
                     c("character","double","double","numeric","numeric","integer","double","integer","integer"))

# PRS 
pdf("/media/je/JE3/R_plot/Plot_PDF_PNG/Savgol_PRS.pdf", width = 8.98, height = 5.86, title = "Résultat Lissage")

foreach(i=1:10) %do%{
  point <- prs[which(prs$Point_id==i),]
  
  p <- ggplot(data=point) + geom_point(aes(x = as.Date(Date,format="%Y%m%d"), y=Original),na.rm = TRUE) +
    geom_line(aes(x = as.Date(Date,format="%Y%m%d"), y=Savgol), color="maroon3") +
    facet_grid (window ~ order)+ theme_bw()+
    scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
    ylab("Valeur")
  
  ggsave(paste("/media/je/JE3/R_plot/Plot_PDF_PNG/savgol/savgol_prs_",i,".png",sep=""))
  
}

dev.off()

#PRS cor
pdf("/media/je/JE3/R_plot/Plot_PDF_PNG/Savgol_PRScor.pdf", width = 8.98, height = 5.86, title = "Résultat Lissage")

foreach(i=1:10) %do%{
  pointcor <- prs_cor[which(prs_cor$Point_id==i),]
  
  p <- ggplot(data=pointcor) + geom_point(aes(x = as.Date(Date,format="%Y%m%d"), y=Original),na.rm = TRUE) +
    geom_line(aes(x = as.Date(Date,format="%Y%m%d"), y=Savgol), color="maroon3") +
    facet_grid (window ~ order)+ theme_bw()+
    scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
    ylab("Valeur")
  
  ggsave(paste("/media/je/JE3/R_plot/Plot_PDF_PNG/savgol/savgol_prscor_",i,".png",sep=""))
  
}

dev.off()