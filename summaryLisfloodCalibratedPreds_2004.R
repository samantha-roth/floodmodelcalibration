#average the predictions
rm(list=ls())

nThin= 500

library(raster)

RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/2004flood/Extent/RunTrue2004_1.asc")
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m<- raster::extract(RunTrue.10m,coords.10m)
rm(RunTrue.10m)


################################################################################
################################################################################
################################################################################

#with no estimate of the discrepancy



################################################################################
#homML no discrepancy

homMLpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)


for(i in 1:nThin){
  #for(i in 458:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/2004flood/homML/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homMLpredMat[,i]<- vals
}

homML_sd<- apply(homMLpredMat,1,sd)
homML_mean<- apply(homMLpredMat,1,mean)
homML_qMat<- t(as.matrix(apply(homMLpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homMLpredMat,homML_sd,homML_mean,homML_qMat,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homML.RData")

beepr::beep()

NAvals<- apply(homMLpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(homMLpredMat[,i])))
}
sum(NAvals)
which(NAvals>0)

rm(homMLpredMat)

################################################################################
#homGP no discrepancy

homGPpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/2004flood/homGP/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPpredMat[,i]<- vals
}

homGP_sd<- apply(homGPpredMat,1,sd)
homGP_mean<- apply(homGPpredMat,1,mean)
homGP_qMat<- t(as.matrix(apply(homGPpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPpredMat, homGP_sd,homGP_mean,homGP_qMat,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homGP.RData")


beepr::beep()

NAvals<- apply(homGPpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(homGPpredMat[,i])))
}
sum(NAvals)
which(NAvals>0)

rm(homGPpredMat)

################################################################################
#homGP no discrepancy

homGPCheappredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/2004flood/homGPCheap/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPCheappredMat[,i]<- vals
}

homGPCheap_sd<- apply(homGPCheappredMat,1,sd)
homGPCheap_mean<- apply(homGPCheappredMat,1,mean)
homGPCheap_qMat<- t(as.matrix(apply(homGPCheappredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPCheappredMat,homGPCheap_sd,homGPCheap_mean,homGPCheap_qMat,
     file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/2004flood/homGPCheap.RData")


beepr::beep()

NAvals<- apply(homGPCheappredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(homGPCheappredMat[,i])))
}
sum(NAvals)
which(NAvals>0)

rm(homGPCheappredMat)
################################################################################

################################################################################