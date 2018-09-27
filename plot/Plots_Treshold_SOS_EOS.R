library (readr)
library(ggplot2)
library(foreach)
library(scales)
library(Metrics)
library(ggthemes)

metrics_hants <- read.csv("/media/je/SAUVEGARDE/Cours_SIGMA/1001_Stage/Stage2018/Process/METRICS/METRICS_CheckTreshold_NDVI_hants_PRScor.csv", colClasses = 
                            c("character","character","character","character","character","character","character","character","character",
                              "character","character","character","character"))

metrics_whit <- read.csv("/media/je/SAUVEGARDE/Cours_SIGMA/1001_Stage/Stage2018/Process/METRICS/METRICS_CheckTreshold_NDVI_whittaker_PRScor.csv", colClasses = 
                           c("character","character","character","character","character","character","character","character","character",
                             "character","character","character","character"))

# Reformat Date columns 

metrics_hants$SOS_10 <- as.Date(metrics_hants$SOS_10,format="%Y%m%d")
metrics_hants$SOS_20 <- as.Date(metrics_hants$SOS_20,format="%Y%m%d")
metrics_hants$SOS_30 <- as.Date(metrics_hants$SOS_30,format="%Y%m%d")
metrics_hants$SOS_50 <- as.Date(metrics_hants$SOS_50,format="%Y%m%d")
metrics_hants$Semi <- as.Date(metrics_hants$Semi,format="%d/%m/%Y")
metrics_hants$EOS_50 <- as.Date(metrics_hants$EOS_50,format="%Y%m%d")
metrics_hants$EOS_60 <- as.Date(metrics_hants$EOS_60,format="%Y%m%d")
metrics_hants$EOS_70 <- as.Date(metrics_hants$EOS_70,format="%Y%m%d")
metrics_hants$EOS_80 <- as.Date(metrics_hants$EOS_80,format="%Y%m%d")
metrics_hants$Recolte <- as.Date(metrics_hants$Recolte,format="%d/%m/%Y")
metrics_hants$Methode <- "HANTS"

metrics_whit$SOS_10 <- as.Date(metrics_whit$SOS_10,format="%Y%m%d")
metrics_whit$SOS_20 <- as.Date(metrics_whit$SOS_20,format="%Y%m%d")
metrics_whit$SOS_30 <- as.Date(metrics_whit$SOS_30,format="%Y%m%d")
metrics_whit$SOS_50 <- as.Date(metrics_whit$SOS_50,format="%Y%m%d")
metrics_whit$Semi <- as.Date(metrics_whit$Semi,format="%d/%m/%Y")
metrics_whit$EOS_50 <- as.Date(metrics_whit$EOS_50,format="%Y%m%d")
metrics_whit$EOS_60 <- as.Date(metrics_whit$EOS_60,format="%Y%m%d")
metrics_whit$EOS_70 <- as.Date(metrics_whit$EOS_70,format="%Y%m%d")
metrics_whit$EOS_80 <- as.Date(metrics_whit$EOS_80,format="%Y%m%d")
metrics_whit$Recolte <- as.Date(metrics_whit$Recolte,format="%d/%m/%Y")
metrics_whit$Methode <- "Whittaker"

# Bind data
metrics <- rbind(metrics_hants,metrics_whit)
# metrics <- subset(metrics, metrics$Projet=="Oracle")

metrics$SystCult[metrics$SystCult=="Mixed_Groundnut"] <- "Arachide Mixte"
metrics$SystCult[metrics$SystCult=="Mixed_Millet"] <- "Mil Mixte"
metrics$SystCult[metrics$SystCult=="Pure_Millet"] <- "Mil Pur"

metrics$CropSyst2[metrics$SystCult=="Arachide Mixte"] <- paste("Arachide Mixte ","(",
                                                               length(which(metrics$SystCult=="Arachide Mixte"))/2, ")", sep="")
metrics$CropSyst2[metrics$SystCult=="Mil Mixte"] <- paste("Mil Mixte ","(",
                                                          length(which(metrics$SystCult=="Mil Mixte"))/2, ")", sep="")
metrics$CropSyst2[metrics$SystCult=="Mil Pur"] <- paste("Mil Pur ","(",
                                                        length(which(metrics$SystCult=="Mil Pur"))/2, ")", sep="")

# Get each Crop Syst data

treshold <- c("10","20","30","50")
methode <- c("HANTS","Whittaker")
SystCult <- c("Arachide Mixte","Mil Mixte","Mil Pur")

cv <- function(x){(sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))*100}

# Statistical indicators

indicators <- data.frame(stringsAsFactors = FALSE)

for(i in 1:length(treshold)){
  for (j in 1:length(methode)){
    for (k in 1:length(SystCult)){
      
      # Seuil
      seuil = paste("SOS",treshold[i],sep="_")
      
      # Compute RMSE
      rmse_value <- rmse(as.numeric(format(metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
                         ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      # MAE
      # mae_value <- mae(as.numeric(format(metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
      #                  ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      # Mean
      # mean_value <- mean(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
      #                                     metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))
      
      
      # Std
      std_value <- sd(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
                                          metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))

      # Coefficient of Variation
      cv_value <- cv(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
                                         metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))
      
      # Correlation
      # cor_value <- cor(as.numeric(format(metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
      #                  ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      df <- data.frame(Treshold = treshold[i], Methode = methode[j], stringsAsFactors = FALSE)
      df$CropSyst <- unique(metrics$CropSyst2[metrics$SystCult==SystCult[k]])
      df$RMSE <- rmse_value
      df$Std <- std_value
      df$CV <- cv_value
      # df$Mean <- mean_value
      # df$MAE <- mae_value
      # df$Corr <- cor_value
      
      indicators <- rbind (indicators, df)
      
    }
  }
}

### Plot ###

# RMSE vs Std
p <- ggplot (indicators,aes(x=RMSE, y=Std, color=Methode, shape=Treshold)) + geom_point(size=2) +
  facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
ggsave("C:/Users/User/Desktop/Plots_new/SOS_Oracle_RMSE_vs_Std.png", width=7.18,height=5)

# RMSE vs CV
p <- ggplot (indicators,aes(x=RMSE, y=CV, color=Methode, shape=Treshold)) + geom_point(size=3) + 
  facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw() + 
  xlab("RMSE (jours)") +  ylab("CV (%)") + 
  labs(shape = "Seuil", color = "Méthode")
ggsave("/home/je/Bureau/SOS_RMSE_vs_CV.pdf", width=7.18,height=5)

# # RMSE vs MAE
# p <- ggplot (indicators,aes(x=RMSE, y=MAE, color=Methode, shape=Treshold)) + geom_point(size=2) + 
#   facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
# ggsave("C:/Users/User/Desktop/Plots_new/SOS_Oracle_RMSE_vs_MAE.png", width=7.18,height=5)
# 
# # RMSE vs Cor
# p <- ggplot (indicators,aes(x=RMSE, y=Corr, color=Methode, shape=Treshold)) + geom_point(size=2) + 
#   facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
# ggsave("C:/Users/User/Desktop/Plots_new/SOS_Oracle_RMSE_vs_Corr.png", width=7.18,height=5)


# Retriev Shift between SOS and Sowing Date

diff <- data.frame(stringsAsFactors = FALSE)

for(i in 1:length(treshold)){
  for (j in 1:length(methode)){
    for (k in 1:length(SystCult)){
      
      seuil = paste("SOS",treshold[i],sep="_")
      subtraction <- difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
                              metrics$Semi[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")
      
      df <- data.frame(Ecart=subtraction,stringsAsFactors = FALSE)
      df$Treshold <- treshold[i]
      df$Methode <- methode[j]
      df$CropSyst <- unique(metrics$CropSyst2[metrics$SystCult==SystCult[k]])
      # paste(unlist(strsplit(SystCult[k],"_"))[1],unlist(strsplit(SystCult[k],"_"))[2],sep=" ")
      
      diff <- rbind (diff, df)
    }
  }
}

### BoxPlot ###
p <- ggplot (diff, aes(x=factor(Treshold), y=as.numeric(Ecart), fill=Methode)) + geom_boxplot() + facet_grid(.~ CropSyst) +
  theme_bw() + xlab("Seuils (%)") +  ylab("SOS - Dates de semis (jours)") + 
  labs(fill = "Méthode")
ggsave("/home/je/Bureau/SOS_Boxplot_Oracle.pdf", width=7.18,height=5)


#####################################################################################################################################

# EOS
treshold_eos <- c("50","60","70","80")

# Statistical indicators

indicators_eos <- data.frame(stringsAsFactors = FALSE)

for(i in 1:length(treshold_eos)){
  for (j in 1:length(methode)){
    for (k in 1:length(SystCult)){
      
      # Seuil
      seuil = paste("EOS",treshold_eos[i],sep="_")
      
      # Compute RMSE
      rmse_value <- rmse(as.numeric(format(metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
                         ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      # MAE
      # mae_value <- mae(as.numeric(format(metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
      #                  ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      # Mean
      # mean_value <- mean(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
      #                                        metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))
      
      
      # Std
      # std_value <- sd(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
      #                                     metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))
      
      # Coefficient of Variation
      cv_value <- cv(as.numeric(difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
                                         metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")))
      
      # Correlation
      # cor_value <- cor(as.numeric(format(metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],"%j"))
      #                  ,as.numeric(format(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],"%j")))
      
      df <- data.frame(Treshold = treshold_eos[i], Methode = methode[j], stringsAsFactors = FALSE)
      df$CropSyst <- unique(metrics$CropSyst2[metrics$SystCult==SystCult[k]])
      df$RMSE <- rmse_value
      # df$Std <- std_value
      df$CV <- cv_value
      # df$Mean <- mean_value
      # df$MAE <- mae_value
      # df$Corr <- cor_value
      
      indicators_eos <- rbind (indicators_eos, df)
      
    }
  }
}

### Plot ###

# RMSE vs Std
# p <- ggplot (indicators_eos,aes(x=RMSE, y=Std, color=Methode, shape=Treshold)) + geom_point(size=2) + 
#   facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
# ggsave("C:/Users/User/Desktop/Plots_new/EOS_RMSE_vs_Std.png", width=7.18,height=5)

# RMSE vs CV
p <- ggplot (indicators_eos,aes(x=RMSE, y=CV, color=Methode, shape=Treshold)) + geom_point(size=3) + 
  facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()+ 
  xlab("RMSE (jours)") +  ylab("CV (%)") + labs(shape = "Seuil", color = "Méthode")
ggsave("/home/je/Bureau/EOS_RMSE_vs_CV.pdf", width=7.18,height=5)

# RMSE vs MAE
# p <- ggplot (indicators_eos,aes(x=RMSE, y=MAE, color=Methode, shape=Treshold)) + geom_point(size=2) + 
#   facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
# ggsave("C:/Users/User/Desktop/Plots_new/EOS_RMSE_vs_MAE.png", width=7.18,height=5)

# RMSE vs Cor
# p <- ggplot (indicators_eos,aes(x=RMSE, y=Corr, color=Methode, shape=Treshold)) + geom_point(size=2) + 
#   facet_grid(.~ CropSyst) + scale_shape_manual(values=c(19,17,15,3)) + theme_bw()
# ggsave("C:/Users/User/Desktop/Plots_new/EOS_RMSE_vs_Corr.png", width=7.18,height=5)


# Retriev Shift between SOS and Sowing Date

diff_eos <- data.frame(stringsAsFactors = FALSE)

for(i in 1:length(treshold_eos)){
  for (j in 1:length(methode)){
    for (k in 1:length(SystCult)){
      
      seuil = paste("EOS",treshold_eos[i],sep="_")
      subtraction <- difftime(metrics[seuil][which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k]),],
                              metrics$Recolte[which(metrics$Methode==methode[j] & metrics$SystCult==SystCult[k])],units="days")
      
      df <- data.frame(Ecart=subtraction,stringsAsFactors = FALSE)
      df$Treshold <- treshold_eos[i]
      df$Methode <- methode[j]
      df$CropSyst <- unique(metrics$CropSyst2[metrics$SystCult==SystCult[k]])
      
      diff_eos <- rbind (diff_eos, df)
    }
  }
}

### BoxPlot ###
p <- ggplot (diff_eos, aes(x=factor(Treshold), y=as.numeric(Ecart), fill=Methode)) + geom_boxplot() + facet_grid(.~ CropSyst) +
  theme_bw() + xlab("Seuils (%)") +  ylab("EOS - Dates de récolte (jours)") + 
  labs(fill = "Méthode")
ggsave("/home/je/Bureau/EOS_Boxplot_Oracle.pdf", width=7.18,height=5)
