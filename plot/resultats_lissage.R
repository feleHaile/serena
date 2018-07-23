library(readr)
library(ggplot2)
library(ggthemes)
library(foreach)
library(scales)

# PRS Cor
# NDVI
ndvi <- read.csv("/media/je/JE3/R_plot/Plot_CSV/agg_NDVI_original_PRScor_toggplot.csv",colClasses = 
                   c("Date","character","character","character","character","Date","Date","double","double"))
ndvi <- ndvi[complete.cases(ndvi), ]
ndvi$CropSys[ndvi$CropSys=="Mixed"] <- "Mixte"
ndvi$CropSys[ndvi$CropSys=="Pure"] <- "Pur"
id <- unique(ndvi$ID)

hants <- read.csv("/media/je/JE3/R_plot/Plot_CSV/agg_NDVI_hants_PRScor_toggplot.csv",colClasses = 
                    c("Date","double","double","character"))
hants <- hants[complete.cases(hants),]
hants$Methode <- "HANTS"

whit <- read.csv("//media/je/JE3/R_plot/Plot_CSV/agg_NDVI_whittaker_PRScor_toggplot.csv", colClasses = 
                   c("Date","double","double","character"))
whit <- whit[complete.cases(whit),]
whit$Methode <- "Whittaker"

pdf("/media/je/JE3/R_plot/Plot_PDF_PNG/Lissage_def.pdf", width = 8.98, height = 5.86, title = "Résultat Lissage")

foreach(i=1:length(id)) %do% {
  plot <- ndvi[which(ndvi$ID==id[i]),]
  plot_hants <- hants[which(hants$ID==id[i]),]
  plot_whit <- whit[which(whit$ID==id[i]),]
  lissdf <- rbind(plot_hants,plot_whit)
  
  p <- ggplot() + 
    geom_vline(xintercept = unique(plot$Semi), color="mediumorchid1", size=1.2)+
    geom_vline(xintercept = unique(plot$Recolte), color="mediumorchid1", size=1.2)+
    geom_line(data=lissdf, aes(x=Date,y=Avg, color = Methode), size = 1) + 
    geom_point(data = plot, aes(x= Date, y = Avg, shape = Sensor), size = 2, na.rm = TRUE) +
    ggtitle(paste("Parcelle",plot$ID, "-", plot$Culture, plot$CropSys, sep=" ")) +
    scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
    ylab("Valeur")
  
  p
}

dev.off()

###########"Facets"###########

ndvi_subset <- subset(ndvi,ndvi$ID=="AMES1"|ndvi$ID=="AMSF1"|ndvi$ID=="JA1"|ndvi$ID=="MAEE3"|ndvi$ID=="MMFF1"|ndvi$ID=="P16")
ndvi_subset$Systeme <- paste(ndvi_subset$Culture,ndvi_subset$CropSys, sep=" ")
plot_hants <- subset(hants,hants$ID=="AMES1"|hants$ID=="AMSF1"|hants$ID=="JA1"|hants$ID=="MAEE3"|hants$ID=="MMFF1"|hants$ID=="P16")
plot_whit <- subset(whit,hants$ID=="AMES1"|whit$ID=="AMSF1"|whit$ID=="JA1"|whit$ID=="MAEE3"|whit$ID=="MMFF1"|whit$ID=="P16")

lissdf <- rbind(plot_hants,plot_whit)


p <- ggplot(data=ndvi_subset, aes(x=Date, y=Avg)) + #geom_line(aes(color = Systeme)) + 
  geom_point(aes(shape = Sensor)) +
  geom_vline(aes(xintercept = Semi), colour="grey29", size=1) + #maroon3
  geom_vline(aes(xintercept = Recolte),  colour="grey29", size=1) + #darkorchid3
  #geom_ribbon(aes(ymin=Avg-Std,ymax=Avg+Std, fill=Systeme),alpha=0.2)+
  geom_line(data=lissdf, aes(color = Methode)) +
  facet_wrap (~ ID, ncol = 3) + theme_bw()+
  #ggtitle("Moyenne parcellaire de NDVI et Ecart type") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
  ylab("Valeur") + 
  labs(shape = "Capteurs", color = "Méthode de lissage")

ggsave("/home/je/Bureau/serena/report/materiels_methodes/lissage_prscor.png")
##############################

# MSAVI2

########################################################################################################

# PRS
# NDVI 

ndvi <- read.csv("/media/je/JE3/R_plot/Plot_CSV/agg_NDVI_original_PRS_toggplot.csv",colClasses = 
                   c("Date","character","character","character","character","Date","Date","double","double"))
ndvi <- ndvi[complete.cases(ndvi), ]
ndvi$CropSys[ndvi$CropSys=="Mixed"] <- "Mixte"
ndvi$CropSys[ndvi$CropSys=="Pure"] <- "Pur"
id <- unique(ndvi$ID)

hants <- read.csv("/media/je/JE3/R_plot/Plot_CSV/agg_NDVI_hants_PRS_toggplot.csv",colClasses = 
                    c("Date","double","double","character"))
hants <- hants[complete.cases(hants),]
hants$Methode <- "HANTS"

whit <- read.csv("/media/je/JE3/R_plot/Plot_CSV/agg_NDVI_whittaker_PRS_toggplot.csv", colClasses = 
                   c("Date","double","double","character"))
whit <- whit[complete.cases(whit),]
whit$Methode <- "Whittaker"

pdf("/media/je/JE3/R_plot/Plot_PDF_PNG/Lissage_test2.pdf", width = 8.98, height = 5.86, title = "Résultat Lissage")

foreach(i=1:length(id)) %do% {
  plot <- ndvi[which(ndvi$ID==id[i]),]
  plot_hants <- hants[which(hants$ID==id[i]),]
  plot_whit <- whit[which(whit$ID==id[i]),]
  lissdf <- rbind(plot_hants,plot_whit)
  
  p <- ggplot() + 
    geom_vline(xintercept = unique(plot$Semi), color="mediumorchid1", size=1.2)+
    geom_vline(xintercept = unique(plot$Recolte), color="mediumorchid1", size=1.2)+
    geom_line(data=lissdf, aes(x=Date,y=Avg, color = Methode), size = 1) + 
    geom_point(data = plot, aes(x= Date, y = Avg, shape = Sensor), size = 2, na.rm = TRUE) +
    ggtitle(paste("Parcelle",plot$ID, "-", plot$Culture, plot$CropSys, sep=" ")) +
    scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
    ylab("Valeur")
  
  p
}

dev.off()

###########"Facets"###########

ndvi_subset <- subset(ndvi,ndvi$ID=="AMES1"|ndvi$ID=="AMSF1"|ndvi$ID=="JA1"|ndvi$ID=="MAEE3"|ndvi$ID=="MMFF1"|ndvi$ID=="P16")
ndvi_subset$Systeme <- paste(ndvi_subset$Culture,ndvi_subset$CropSys, sep=" ")
plot_hants <- subset(hants,hants$ID=="AMES1"|hants$ID=="AMSF1"|hants$ID=="JA1"|hants$ID=="MAEE3"|hants$ID=="MMFF1"|hants$ID=="P16")
plot_whit <- subset(whit,hants$ID=="AMES1"|whit$ID=="AMSF1"|whit$ID=="JA1"|whit$ID=="MAEE3"|whit$ID=="MMFF1"|whit$ID=="P16")

lissdf <- rbind(plot_hants,plot_whit)


p <- ggplot(data=ndvi_subset, aes(x=Date, y=Avg)) + #geom_line(aes(color = Systeme)) + 
  geom_point(aes(shape = Sensor)) +
  geom_vline(aes(xintercept = Semi), colour="grey29", size=1) + #maroon3
  geom_vline(aes(xintercept = Recolte),  colour="grey29", size=1) + #darkorchid3
  #geom_ribbon(aes(ymin=Avg-Std,ymax=Avg+Std, fill=Systeme),alpha=0.2)+
  geom_line(data=lissdf, aes(color = Methode)) +
  facet_wrap (~ ID, ncol = 3) + theme_bw() +
  #ggtitle("Moyenne parcellaire de NDVI et Ecart type") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
  ylab("Valeur") + 
  labs(shape = "Capteurs", color = "Méthode de lissage")

ggsave("/home/je/Bureau/serena/report/materiels_methodes/lissage_prs.png")
