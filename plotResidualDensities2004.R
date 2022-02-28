
rm(list=ls())
#plot Lisflood calibration preds
library(raster)
library(RColorBrewer)
library(rgdal)
library(sf)

RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/2004flood/Extent/RunTrue2004_1.asc")
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m<- raster::extract(RunTrue.10m,coords.10m)
########################################################################################################################

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homGP.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homML.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homGPCheap.RData")

########################################################################################################################


Run_homML_mean<- RunTrue.10m
values(Run_homML_mean)<- homML_mean

Run_homGP_mean<- RunTrue.10m
values(Run_homGP_mean)<- homGP_mean

Run_homGPCheap_mean<- RunTrue.10m
values(Run_homGPCheap_mean)<- homGPCheap_mean

#Pdf of all residuals

#Pdf of residuals at the hypothetical houses


#plot(density(homML_mean - truevals.10m), col= "red", main= "PDF of residuals for 2004 flood", xlim= c(-.8, 2.2))
#lines(density(homGP_mean - truevals.10m), col= "blue")
#lines(density(homGPCheap_mean - truevals.10m), col= "turquoise")
#legend(.5,8 , legend=c("HomMR Residual", 
#                         "HomGP10 Residual",
#                         "HomGP50 Residual"),
#       col=c("red",
#             "blue", 
#             "turquoise"), lty=1,cex=1.5,box.lty = 0.1)



residual_df<- data.frame(
  "Residuals"= c(homML_mean- truevals.10m, homGP_mean- truevals.10m, homGPCheap_mean- truevals.10m),
  "Source"= c(rep("HomMR",length(homML_mean)),
              rep("HomGP10",length(homGP_mean)),
              rep("HomGP50",length(homGPCheap_mean)))
)

library(ggplot2)
jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/residualPDF_2004.jpeg",width = 800, height = 600)
ggplot(residual_df, aes(x=Residuals, color=Source))+
  geom_density() + 
  xlab("Residual") + ylab("Density") +
  ggtitle("Residual densities") +
  scale_color_manual(values=c("blue", "turquoise", "red"))+ 
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24))
dev.off()



summary(homML_mean - truevals.10m)
summary(homGP_mean - truevals.10m)
summary(homGPCheap_mean - truevals.10m)


length(which( (homML_mean - truevals.10m) >1))/length(homML_mean - truevals.10m)
length(which( (homGP_mean - truevals.10m) >1))/length(homGP_mean - truevals.10m)
length(which( (homGPCheap_mean - truevals.10m) >1))/length(homGPCheap_mean - truevals.10m)

length(which( (homML_mean - truevals.10m) <0))/length(homML_mean - truevals.10m)
length(which( (homGP_mean - truevals.10m) <0))/length(homGP_mean - truevals.10m)
length(which( (homGPCheap_mean - truevals.10m) <0))/length(homGPCheap_mean - truevals.10m)

