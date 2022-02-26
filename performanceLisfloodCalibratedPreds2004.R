rm(list=ls())
#plot Lisflood calibration preds
library(raster)

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


#par(mfrow=c(3,2))
#plot(RunTrue.10m, main= "Observation")
#plot(Run_hetML_mean, main= "Mean prediction from HetML emulation-calibration")
#plot(Run_homML_mean, main= "Mean prediction from HomML emulation-calibration")
#plot(Run_hetGP_mean, main= "Mean prediction from HetGP emulation-calibration")
#plot(Run_homGP_mean, main= "Mean prediction from HomGP emulation-calibration")


cuts=seq(0,9,by=.01)
pal <- colorRampPalette(c("white","blue"))
plot(RunTrue.10m, breaks=cuts, col = pal(length(cuts)+1)) #plot with defined breaks


plot(RunTrue.10m, main= "Observation", breaks=cuts, col = pal(length(cuts)+1))
par(mfrow=c(3,2))

plot(Run_homML_mean, main= "Mean prediction from HomMR emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
plot(Run_homGP_mean, main= "Mean prediction from HomGP10 emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))
plot(Run_homGPCheap_mean, main= "Mean prediction from HomGP50 emulation-calibration", breaks=cuts, col = pal(length(cuts)+1))


Run_homML_diff<- RunTrue.10m
values(Run_homML_diff)<- homML_mean - truevals.10m

Run_homGP_diff<- RunTrue.10m
values(Run_homGP_diff)<- homGP_mean - truevals.10m

Run_homGPCheap_diff<- RunTrue.10m
values(Run_homGPCheap_diff)<- homGPCheap_mean - truevals.10m


maxdiff<- max(c(homML_mean - truevals.10m, homGP_mean - truevals.10m, homGPCheap_mean - truevals.10m))
mindiff<- min(c(homML_mean - truevals.10m, homGP_mean - truevals.10m, homGPCheap_mean - truevals.10m))


cutsdiff=seq(round(mindiff,1),round(maxdiff,1),by=.01)
pal <- colorRampPalette(c("white","blue"))

#plot(RunTrue.10m)
par(mfrow=c(2,2))
plot(RunTrue.10m, main= "Observation", breaks=cuts, col = pal(length(cuts)+1))
plot(Run_homML_diff, main= "Mean prediction from HomML emulation-calibration minus observation", breaks=cutsdiff, col = pal(length(cutsdiff)+1))
plot(Run_homGP_diff, main= "Mean prediction from HomGP10 emulation-calibration minus observation", breaks=cutsdiff, col = pal(length(cutsdiff)+1))
plot(Run_homGPCheap_diff, main= "Mean prediction from HomGP50 emulation-calibration minus observation", breaks=cutsdiff, col = pal(length(cutsdiff)+1))



Run_homML_sd<- RunTrue.10m
values(Run_homML_sd)<- homML_sd


Run_homGP_sd<- RunTrue.10m
values(Run_homGP_sd)<- homGP_sd

Run_homGPCheap_sd<- RunTrue.10m
values(Run_homGPCheap_sd)<- homGPCheap_sd

par(mfrow=c(2,2))
plot(Run_homML_sd, main= "HomML sd")
plot(Run_homGP_sd, main= "HomGP sd")
plot(Run_homGPCheap_sd, main= "HomGPCheap sd")
#see less extreme range of variation in predictions from HetML and HomML 
#compared to HetGP
#HomGP looks similar


################################################################################
#Calculate Euclidean Distances

homML_ED<- sqrt(sum((homML_mean - truevals.10m)^2))
homGP_ED<- sqrt(sum((homGP_mean - truevals.10m)^2))
homGPCheap_ED<- sqrt(sum((homGPCheap_mean - truevals.10m)^2))

homML_ED; homGP_ED; homGPCheap_ED

#RMSE

homML_RMSE<- sqrt(mean((homML_mean - truevals.10m)^2))
homGP_RMSE<- sqrt(mean((homGP_mean - truevals.10m)^2))
homGPCheap_RMSE<- sqrt(mean((homGPCheap_mean - truevals.10m)^2))

homML_RMSE; homGP_RMSE; homGPCheap_RMSE

#PBIAS


homML_PBias<- 100 * (sum( homML_mean - truevals.10m ) / sum( truevals.10m ))  
homGP_PBias<- 100 * (sum( homGP_mean - truevals.10m ) / sum( truevals.10m ))  
homGPCheap_PBias<- 100 * (sum( homGPCheap_mean - truevals.10m ) / sum( truevals.10m ))  


homML_PBias
homGP_PBias 
homGPCheap_PBias 

#smaller Euclidean Distances from HetML and HomML than HetGP and HomGP

################################################################################

#Calculate Fit and Correctness
res.e<- 10

mGridArea<- res.e^2
rGridArea<- res.e^2

Ar<- length(which(truevals.10m>0))*rGridArea


Am<- length(which(homML_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homML_mean>0)))*rGridArea
Fvals_homML<- Arm/(Am + Ar - Arm)
Cvals_homML<- Arm/Ar

Am<- length(which(homGP_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homGP_mean>0)))*rGridArea
Fvals_homGP<- Arm/(Am + Ar - Arm)
Cvals_homGP<- Arm/Ar

Am<- length(which(homGPCheap_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homGPCheap_mean>0)))*rGridArea
Fvals_homGPCheap<- Arm/(Am + Ar - Arm)
Cvals_homGPCheap<- Arm/Ar


Fvals_homML
Fvals_homGP
Fvals_homGPCheap


Cvals_homML
Cvals_homGP
Cvals_homGPCheap


save(homML_ED, homGP_ED, homGPCheap_ED, Fvals_homML, Fvals_homGP, Fvals_homGPCheap, Cvals_homML, Cvals_homGP, Cvals_homGPCheap,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/ED_F_Cvals.RData")


#smaller Euclidean Distances from HetML and HomML than HetGP and HomGP
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/ED_F_Cvals_delta.RData")
################################################################################


#look at results from approaches using only cheap model runs


Run_hetGPCheapD_diff<- RunTrue.10m
values(Run_hetGPCheapD_diff)<- hetGPCheapD_mean - truevals.10m

Run_homGPCheapD_diff<- RunTrue.10m
values(Run_homGPCheapD_diff)<- homGPCheapD_mean - truevals.10m


#plot(RunTrue.10m)
par(mfrow=c(1,2))
plot(Run_hetGPCheapD_diff, main= "Mean prediction from HetGP (Cheap) emulation-calibration minus observation")
plot(Run_homGPCheapD_diff, main= "Mean prediction from HomGP (Cheap) emulation-calibration minus observation")

#sd
Run_hetGPCheapD_sd<- RunTrue.10m
values(Run_hetGPCheapD_sd)<- hetGPCheapD_sd

Run_homGPCheapD_sd<- RunTrue.10m
values(Run_homGPCheapD_sd)<- homGPCheapD_sd

par(mfrow=c(1,2))
plot(Run_hetGPCheapD_sd, main= "HetGPCheapD sd")
plot(Run_homGPCheapD_sd, main= "HomGPCheapD sd")
#see less extreme range of variation in predictions from HetML and HomML 
#compared to HetGP
#HomGP looks similar


################################################################################
#Calculate Euclidean Distances

hetGPCheap_ED<- sqrt(sum((hetGPCheap_mean - truevals.10m)^2))
homGPCheap_ED<- sqrt(sum((homGPCheap_mean - truevals.10m)^2))

hetGPCheap_ED; homGPCheap_ED

hetGPCheap_RMSE<- sqrt(mean((hetGPCheap_mean - truevals.10m)^2))
homGPCheap_RMSE<- sqrt(mean((homGPCheap_mean - truevals.10m)^2))

hetGPCheap_RMSE; homGPCheap_RMSE
#smaller Euclidean Distances from HetML and HomML than HetGP and HomGP

hetGPCheap_PBias<- 100 * (sum( hetGPCheap_mean - truevals.10m ) / sum( truevals.10m ))  
homGPCheap_PBias<- 100 * (sum( homGPCheap_mean - truevals.10m ) / sum( truevals.10m ))  

hetGPCheap_PBias
homGPCheap_PBias 

################################################################################

#Calculate Fit and Correctness
res.e<- 10

mGridArea<- res.e^2
rGridArea<- res.e^2

Ar<- length(which(truevals.10m>0))*rGridArea

Am<- length(which(hetGPCheap_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(hetGPCheap_mean>0)))*rGridArea
Fvals_hetGPCheap<- Arm/(Am + Ar - Arm)
Cvals_hetGPCheap<- Arm/Ar

Am<- length(which(homGPCheap_mean>0))*mGridArea
Arm<- length(intersect(which(truevals.10m>0),which(homGPCheap_mean>0)))*rGridArea
Fvals_homGPCheap<- Arm/(Am + Ar - Arm)
Cvals_homGPCheap<- Arm/Ar


Fvals_hetGPCheap
Fvals_homGPCheap

Cvals_hetGPCheap
Cvals_homGPCheap

save(hetGPCheap_ED, homGPCheap_ED, 
     Fvals_hetGPCheap, Fvals_homGPCheap, 
     Cvals_hetGPCheap, Cvals_homGPCheap, 
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/ED_F_Cvals_delta_cheap.RData")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/ED_F_Cvals_delta_cheap.RData")


################################################################################
#look at these metrics for the median


Run_hetML_mean<- RunTrue.10m
values(Run_hetML_mean)<- hetML_qMat[,3]

Run_homML_mean<- RunTrue.10m
values(Run_homML_mean)<- homML_qMat[,3]

Run_homGP_mean<- RunTrue.10m
values(Run_homGP_mean)<- homGP_qMat[,3]

