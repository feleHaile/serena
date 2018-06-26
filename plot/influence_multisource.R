library(readr)
library(ggplot2)
library(ggthemes)
library(foreach)
library(scales)

ndvi <- read.csv("/media/je/JE3/to_plot/Plot_CSV/agg_NDVI_original_PRS_toggplot.csv",colClasses = 
                   c("Date","character","character","character","character","Date","Date","double","double"))
ndvi$CropSys[ndvi$CropSys=="Mixed"] <- "Mixte"
ndvi <- ndvi[complete.cases(ndvi), ]
ndvi$index <- "NDVI"

msavi2 <- read.csv("/media/je/JE3/to_plot/Plot_CSV/agg_MSAVI2_original_PRS_toggplot.csv",colClasses = 
                     c("Date","character","character","character","character","Date","Date","double","double"))
msavi2$CropSys[msavi2$CropSys=="Mixed"] <- "Mixte"
msavi2 <- msavi2[complete.cases(msavi2), ]
msavi2$index <- "MSAVI2"

id <- unique(ndvi$ID)

ndvi_subset <- subset(ndvi,ndvi$ID=="MAEE3"|ndvi$ID=="MAEE2"|ndvi$ID=="P16"|ndvi$ID=="P12")
ndvi_subset$Systeme <- paste(ndvi_subset$Culture,ndvi_subset$CropSys, sep=" ")

###########"Facets"###########
p <- ggplot(data=ndvi_subset, aes(x=Date, y=Avg)) + geom_point(aes(shape = Sensor, color = Systeme)) + geom_line(aes(color = Systeme)) + 
  geom_vline(xintercept = unique(plot$Semi), color="darkorange", size=1.2)+geom_vline(xintercept = unique(plot$Recolte), color="mediumorchid1", size=1.2) +
  facet_wrap (~ ID, ncol = 3)

p
##############################

pdf("/media/je/JE3/to_plot/Plot_PDF_PNG/test.pdf", width = 8.98, height = 5.86, title = "Moyenne Parcellaire NDVI")

foreach(i=1:length(id)) %do% {
  plot <- ndvi[which(ndvi$ID==id[i]),]
  p <- ggplot(data=plot, aes(x=Date, y=Avg)) + 
    geom_vline(xintercept = unique(plot$Semi), color="darkorange", size=1.2)+geom_vline(xintercept = unique(plot$Recolte), color="darkorange", size=1.2) +
    geom_point(aes(shape = Sensor)) +
    geom_line() +
    geom_ribbon(aes(ymin=plot$Avg-plot$Std,ymax=plot$Avg+plot$Std),alpha=0.2) + 
    ggtitle(paste("Parcelle",plot$ID, "-", plot$Culture, plot$CropSys, sep=" ")) +
    xlab("Temps") + ylab("Valeur")

  print (p)
}
dev.off()

pdf("/media/je/JE3/to_plot/Plot_PDF_PNG/test2.pdf", width = 8.98, height = 5.86, title = "Moyenne Parcellaire MSAVI2")

foreach(i=1:length(id)) %do% {
  plot <- msavi2[which(msavi2$ID==id[i]),]
  p <- ggplot(data=plot, aes(x=Date, y=Avg)) + 
    geom_vline(xintercept = unique(plot$Semi), color="darkorange", size=1.2)+geom_vline(xintercept = unique(plot$Recolte), color="darkorange", size=1.2) +
    geom_point(aes(shape = Sensor)) +
    geom_line() +
    geom_ribbon(aes(ymin=plot$Avg-plot$Std,ymax=plot$Avg+plot$Std),alpha=0.2) + 
    ggtitle(paste("Parcelle",plot$ID, "-", plot$Culture, plot$CropSys, sep=" ")) +
    scale_x_date(name="Date", breaks = date_breaks("2 months"),labels = date_format("%b"), limits = c(as.Date("2017-05-08"),as.Date("2017-11-19"))) +
    ylab("Valeur")
  
  print (p)
}
dev.off()

pdf("/media/je/JE3/to_plot/Plot_PDF_PNG/test.pdf", width = 8.98, height = 5.86, title = "Moyenne Parcellaire NDVI et MSAVI2")

foreach(i=1:length(id)) %do% {
  plot1 <- ndvi[which(ndvi$ID==id[i]),]
  plot2 <- msavi2[which(ndvi$ID==id[i]),]
  indexdf <- rbind(plot1,plot2)
  
  p <- ggplot() +
    geom_vline(xintercept = unique(indexdf$Semi), color="mediumorchid1", size=1.2)+ 
    geom_vline(xintercept = unique(indexdf$Recolte), color="mediumorchid1", size=1.2) +
    geom_ribbon(data = indexdf, aes(x = Date, ymin=indexdf$Avg-indexdf$Std,ymax=indexdf$Avg+indexdf$Std, fill = index),alpha=0.2)+geom_line(data = indexdf, aes(x=Date, y=Avg, color = index)) + 
    geom_point(data = indexdf, aes(x=Date, y=Avg, shape = Sensor)) + 
    ggtitle(paste("Parcelle",plot$ID, "-", plot$Culture, plot$CropSys, sep=" ")) +
    xlab("Temps") + ylab("Valeur")
  
  print (p)
}
dev.off()