#Joint density plots


#get emulator predictions from MCMC samples
rm(list=ls())
source("C:/FloodingModelCalibrationProject/sml-athena-main/GPfunctionsOptim.R")
source("C:/FloodingModelCalibrationProject/sml-athena-main/hetGPfunctions.R")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP.RData")
res_homGP<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP800.RData")
res_homGP800<- res[-1,]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
res_homML<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
res_homGPCheap<- res[-c(1,(1e5)+1),]



rm(res)

library(ggplot2)
res_homMLdf<- as.data.frame(res_homML)
res_homGPdf<- as.data.frame(res_homGP)
res_homGPCheapdf<- as.data.frame(res_homGPCheap)

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomMR_jointChREE.jpeg",width = 800, height = 600)
ggplot(res_homMLdf, aes(x=x1.out, y=x4.out) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("HomMR- Joint Histogram") + 
  xlab("Channel Roughness") + ylab("Riverbed Elevation Error") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP10_jointChREE.jpeg",width = 800, height = 600)
ggplot(res_homGPdf, aes(x=x1.out, y=x4.out) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("HomGP10- Joint Histogram") +
  xlab("Channel Roughness") + ylab("Riverbed Elevation Error") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/HomGP50_jointChREE.jpeg",width = 800, height = 600)
ggplot(res_homGPCheapdf, aes(x=x1.out, y=x4.out) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("HomGP50- Joint Histogram") + 
  xlab("Channel Roughness") + ylab("Riverbed Elevation Error") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.text= element_text(size=24),
        legend.title= element_text(size=24))
dev.off()



#plot(res_homGP800[,1],res_homGP800[,4])

library(LaplacesDemon)

##########CHANNEL ROUGHNESS AND RIVERBED ELEVATION ERROR

Channel_Roughness<- res_homGP800[,1]
Riverbed_Elevation_Error<- res_homGP800[,4]

joint.density.plot(Channel_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

Channel_Roughness<- res_homGP[,1]
Riverbed_Elevation_Error<- res_homGP[,4]

joint.density.plot(Channel_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

Channel_Roughness<- res_homML[,1]
Riverbed_Elevation_Error<- res_homML[,4]

joint.density.plot(Channel_Roughness, Riverbed_Elevation_Error, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


Channel_Roughness<- res_homGPCheap[,1]
Riverbed_Elevation_Error<- res_homGPCheap[,4]

joint.density.plot(Channel_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")

##########CHANNEL ROUGHNESS AND RIVER WIDTH ERROR

Channel_Roughness<- res_homGP800[,1]
River_Width_Error<- res_homGP800[,3]

joint.density.plot(Channel_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

Channel_Roughness<- res_homGP[,1]
River_Width_Error<- res_homGP[,3]

joint.density.plot(Channel_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

Channel_Roughness<- res_homML[,1]
River_Width_Error<- res_homML[,3]

joint.density.plot(Channel_Roughness, River_Width_Error, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


Channel_Roughness<- res_homGPCheap[,1]
River_Width_Error<- res_homGPCheap[,3]

joint.density.plot(Channel_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")


##########CHANNEL ROUGHNESS AND FLOODPLAIN ROUGHNESS

Channel_Roughness<- res_homGP800[,1]
Floodplain_Roughness<- res_homGP800[,2]

joint.density.plot(Channel_Roughness, Floodplain_Roughness,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

Channel_Roughness<- res_homGP[,1]
Floodplain_Roughness<- res_homGP[,2]

joint.density.plot(Channel_Roughness,Floodplain_Roughness,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

Channel_Roughness<- res_homML[,1]
Floodplain_Roughness<- res_homML[,2]

joint.density.plot(Channel_Roughness, Floodplain_Roughness, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


Channel_Roughness<- res_homGPCheap[,1]
Floodplain_Roughness<- res_homGPCheap[,2]

joint.density.plot(Channel_Roughness, Floodplain_Roughness,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")

##########RIVER WIDTH ERROR AND RIVERBED ELEVATION ERROR

River_Width_Error<- res_homGP800[,3]
Riverbed_Elevation_Error<- res_homGP800[,4]

joint.density.plot(River_Width_Error, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

River_Width_Error<- res_homGP[,3]
Riverbed_Elevation_Error<- res_homGP[,4]

joint.density.plot(River_Width_Error, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

River_Width_Error<- res_homML[,3]
Riverbed_Elevation_Error<- res_homML[,4]

joint.density.plot(River_Width_Error, Riverbed_Elevation_Error, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


River_Width_Error<- res_homGPCheap[,3]
Riverbed_Elevation_Error<- res_homGPCheap[,4]

joint.density.plot(River_Width_Error, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")

##########FLOODPLAIN ROUGHNESS AND RIVERBED ELEVATION ERROR

Floodplain_Roughness<- res_homGP800[,2]
Riverbed_Elevation_Error<- res_homGP800[,4]

joint.density.plot(Floodplain_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

Floodplain_Roughness<- res_homGP[,2]
Riverbed_Elevation_Error<- res_homGP[,4]

joint.density.plot(Floodplain_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

Floodplain_Roughness<- res_homML[,2]
Riverbed_Elevation_Error<- res_homML[,4]

joint.density.plot(Floodplain_Roughness, Riverbed_Elevation_Error, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


Floodplain_Roughness<- res_homGPCheap[,2]
Riverbed_Elevation_Error<- res_homGPCheap[,4]

joint.density.plot(Floodplain_Roughness, Riverbed_Elevation_Error,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")

##########FLOODPLAIN ROUGHNESS AND RIVER WIDTH ERROR

Floodplain_Roughness<- res_homGP800[,2]
River_Width_Error<- res_homGP800[,3]

joint.density.plot(Floodplain_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate (800)")

Floodplain_Roughness<- res_homGP[,2]
River_Width_Error<- res_homGP[,3]

joint.density.plot(Floodplain_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP10- Joint posterior density estimate")

Floodplain_Roughness<- res_homML[,2]
River_Width_Error<- res_homML[,3]

joint.density.plot(Floodplain_Roughness, River_Width_Error, 
                   color= TRUE, 
                   Title= "HomML- Joint posterior density estimate")


Floodplain_Roughness<- res_homGPCheap[,2]
River_Width_Error<- res_homGPCheap[,3]

joint.density.plot(Floodplain_Roughness, River_Width_Error,
                   color= TRUE, 
                   Title= "HomGP50- Joint posterior density estimate")
