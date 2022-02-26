#plot individual parameter density plots


#get emulator predictions from MCMC samples
rm(list=ls())
source("C:/FloodingModelCalibrationProject/sml-athena-main/GPfunctionsOptim.R")
source("C:/FloodingModelCalibrationProject/sml-athena-main/hetGPfunctions.R")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP.RData")
res_homGP<- res[-c(1,(1e5)+1),]

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP800.RData")
#res_homGP800<- res[-1,]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
res_homML<- res[-c(1,(1e5)+1),]

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetSML.RData")
#res_hetML<- res[-c(1,(1e5)+1,(2e5)+1 ),]

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP.RData")
#res_hetGP<- res[-c(1,(1e5)+1,(2e5)+1 ),]

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP_cheap.RData")
#res_hetGPCheap<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
res_homGPCheap<- res[-c(1,(1e5)+1),]

rm(res)

################################################################################
library(ggplot2)

MCMClength<- nrow(res_homML)


rep_homML<- rep("HomMR",MCMClength)
rep_homGP10<- rep("HomGP10",MCMClength)
rep_homGP50<- rep("HomGP50",MCMClength)

allMods<- c(rep_homML, rep_homGP10, rep_homGP50)

#allCh<- cbind(c(res_homML[,1],res_homGP[,1],res_homGPCheap[,1]),allMods)
#allFp<- cbind(c(res_homML[,2],res_homGP[,2],res_homGPCheap[,2]),allMods)
#allRWE<- cbind(c(res_homML[,3],res_homGP[,3],res_homGPCheap[,3]),allMods)
#allREE<- cbind(c(res_homML[,4],res_homGP[,4],res_homGPCheap[,4]),allMods)


allCh<- data.frame(
  "Ch"= c(res_homML[,1],res_homGP[,1],res_homGPCheap[,1]),
  "Source"= allMods
)

allFp<- data.frame(
  "Fp"= c(res_homML[,2],res_homGP[,2],res_homGPCheap[,2]),
  "Source"= allMods
)

allRWE<- data.frame(
  "RWE"= c(res_homML[,3],res_homGP[,3],res_homGPCheap[,3]),
  "Source"= allMods
)


allREE<- data.frame(
  "REE"= c(res_homML[,4],res_homGP[,4],res_homGPCheap[,4]),
  "Source"= allMods
)

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/ChHomPosteriors.jpeg",width = 800, height = 500)
ggplot(allCh, aes(x=Ch, color=Source)) +
  geom_density() + 
  geom_vline(xintercept = .0305, color="black") + 
  xlab("Channel roughness") + ylab("Density") +
  ggtitle("Channel roughness posterior densities") +
  scale_color_manual(values=c("blue", "turquoise", "red"))+ 
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24)) + 
  annotate(geom="text", x=.04, y=1, label="True value",
           color="black", size = 8) +
  geom_segment(x = .02, y = 1/(.1-.02), xend = .1, yend = 1/(.1-.02), color= "dark gray") 
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/FpHomPosteriors.jpeg",width = 800, height = 500)
ggplot(allFp, aes(x=Fp, color=Source))+
  geom_density() + 
  geom_vline(xintercept = .045, color="black") + 
  #geom_hline(yintercept =  1/(.4-.02), color="black") + 
  xlab("Floodplain roughness") + ylab("Density") +
  ggtitle("Floodplain roughness posterior densities") +
  scale_color_manual(values=c("blue", "turquoise", "red"))+ 
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24)) +
  annotate(geom="text", x=.09, y=2, label="True value",
           color="black", size = 8) +
  geom_segment(x = .02, y = 1/(.4-.02), xend = .4, yend = 1/(.4-.02), color= "dark gray")
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/RWEHomPosteriors.jpeg",width = 800, height = 500)
ggplot(allRWE, aes(x=RWE, color=Source)) +
  geom_density() + 
  geom_vline(xintercept = 1, color="black") + 
  #geom_hline(yintercept =  1/(1.05-.95), color="black") + 
  xlab("River width error") + ylab("Density") +
  ggtitle("River width error posterior densities") +
  scale_color_manual(values=c("blue", "turquoise", "red"))+ 
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))+
  annotate(geom="text", x= 1.01, y=2, label="True value",
           color="black", size = 8) +
  geom_segment(x = .95, y = 1/(1.05-.95), xend = 1.05, yend = 1/(1.05-.95), color= "dark gray")
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/calibrationResults/REEHomPosteriors.jpeg",width = 800, height = 500)
ggplot(allREE, aes(x=REE, color=Source)) +
  geom_density() + 
  geom_vline(xintercept = 0, color="black") + 
  #geom_hline(yintercept =  1/(5+5), color="black") + 
  xlab("Riverbed elevation error") + ylab("Density") +
  ggtitle("Riverbed elevation error posterior densities") +
  scale_color_manual(values=c("blue", "turquoise", "red"))+ 
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        legend.text= element_text(size=24))+
  annotate(geom="text", x= -1, y=.05, label="True value",
           color="black", size = 8) +
  geom_segment(x = -5, y = 1/(5+5), xend = 5, yend = 1/(5+5), color= "dark gray")
dev.off()
