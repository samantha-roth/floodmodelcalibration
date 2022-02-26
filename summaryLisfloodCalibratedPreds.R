#average the predictions
rm(list=ls())

nThin= 500

library(raster)

RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m<- raster::extract(RunTrue.10m,coords.10m)
rm(RunTrue.10m)

#with an estimate of delta 

#sample.inds<-seq(10,490,by=20)

#for(i in 400:410){
#  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetGPD/Extent/Run_",i,".asc",sep=""))
#  plot(run)
#} #402 is a problem (but I fixed it)

################################################################################
#hetMLD - with discrepancy

hetMLDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetMLD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetMLDpredMat[,i]<- vals
}

hetMLD_sd<- apply(hetMLDpredMat,1,sd)
hetMLD_mean<- apply(hetMLDpredMat,1,mean)
hetMLD_qMat<- t(as.matrix(apply(hetMLDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(hetMLD_sd,hetMLD_mean,hetMLD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetMLD.RData")

rm(hetMLDpredMat)

################################################################################
#homMLD- with discrepancy

homMLDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homMLD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homMLDpredMat[,i]<- vals
}

homMLD_sd<- apply(homMLDpredMat,1,sd)
homMLD_mean<- apply(homMLDpredMat,1,mean)
homMLD_qMat<- t(as.matrix(apply(homMLDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homMLD_sd,homMLD_mean,homMLD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homMLD.RData")

rm(homMLDpredMat)

################################################################################
#hetGPD- with discrepancy

hetGPDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 402:nThin){
  if(i%%10==0) print(i)
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetGPD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetGPDpredMat[,i]<- vals
}

#NAvals<- apply(hetGPDpredMat,2,function(x) which(is.na(x)))

#NAvals<- rep(NA, nThin)
#for(i in 1:nThin){
#  NAvals[i]<- length(which(is.na(hetGPDpredMat[,i])))
#}

hetGPD_sd<- apply(hetGPDpredMat,1,sd)
hetGPD_mean<- apply(hetGPDpredMat,1,mean)
hetGPD_qMat<- t(as.matrix(apply(hetGPDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))

save(hetGPD_sd,hetGPD_mean,hetGPD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPD.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPD.RData")
rm(homGPDpredMat)

################################################################################
#homGPD- with discrepancy

homGPDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homGPD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPDpredMat[,i]<- vals
}

homGPD_sd<- apply(homGPDpredMat,1,sd)
homGPD_mean<- apply(homGPDpredMat,1,mean)
homGPD_qMat<- t(as.matrix(apply(homGPDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPD_sd,homGPD_mean,homGPD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGPD.RData")

NAvals<- apply(homGPDpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(homGPDpredMat[,i])))
}
sum(NAvals)

rm(homGPDpredMat)

################################################################################

#hetGPCheapD- with discrepancy

hetGPCheapDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetGPCheapD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetGPCheapDpredMat[,i]<- vals
}

hetGPCheapD_sd<- apply(hetGPCheapDpredMat,1,sd)
hetGPCheapD_mean<- apply(hetGPCheapDpredMat,1,mean)
hetGPCheapD_qMat<- t(as.matrix(apply(hetGPCheapDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(hetGPCheapD_sd,hetGPCheapD_mean,hetGPCheapD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPCheapD.RData")

NAvals<- apply(hetGPCheapDpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(hetGPCheapDpredMat[,i])))
}
sum(NAvals)

rm(hetGPCheapDpredMat)


################################################################################

#homGPCheapD- with discrepancy

homGPCheapDpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homGPCheapD/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPCheapDpredMat[,i]<- vals
}

homGPCheapD_sd<- apply(homGPCheapDpredMat,1,sd)
homGPCheapD_mean<- apply(homGPCheapDpredMat,1,mean)
homGPCheapD_qMat<- t(as.matrix(apply(homGPCheapDpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPCheapD_sd,homGPCheapD_mean,homGPCheapD_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGPCheapD.RData")

NAvals<- apply(homGPCheapDpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(homGPCheapDpredMat[,i])))
}
sum(NAvals)

rm(homGPCheapDpredMat)

################################################################################

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPD.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGPD.RData")


identical(hetGPD_mean,homGPD_mean)
identical_qs<- rep(NA, ncol(hetGPD_qMat))
for(i in 1:ncol(hetGPD_qMat)){
  identical_qs[i]<- identical(hetGPD_qMat[,i],homGPD_qMat[,i])
}


################################################################################
################################################################################
################################################################################

#with no estimate of the discrepancy

################################################################################
#hetML no discrepancy
nThin= 500

hetMLpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetML/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetMLpredMat[,i]<- vals
}

hetML_sd<- apply(hetMLpredMat,1,sd)
hetML_mean<- apply(hetMLpredMat,1,mean)
hetML_qMat<- t(as.matrix(apply(hetMLpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(hetMLpredMat,hetML_sd,hetML_mean,hetML_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetML.RData")


NAvals<- apply(hetMLpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(hetMLpredMat[,i])))
}
sum(NAvals)

rm(hetMLpredMat)

################################################################################
#homML no discrepancy

homMLpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)


for(i in 1:nThin){
#for(i in 458:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homML/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homMLpredMat[,i]<- vals
}

homML_sd<- apply(homMLpredMat,1,sd)
homML_mean<- apply(homMLpredMat,1,mean)
homML_qMat<- t(as.matrix(apply(homMLpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homMLpredMat,homML_sd,homML_mean,homML_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homML.RData")

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
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homGP/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPpredMat[,i]<- vals
}

homGP_sd<- apply(homGPpredMat,1,sd)
homGP_mean<- apply(homGPpredMat,1,mean)
homGP_qMat<- t(as.matrix(apply(homGPpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPpredMat, homGP_sd,homGP_mean,homGP_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGP.RData")


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

#hetGP no discrepancy

hetGPpredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){ #147 = problem
#for(i in 148:nThin){
  if((i%%100)==0) print(i)
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetGP/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetGPpredMat[,i]<- vals
}

hetGP_sd<- apply(hetGPpredMat,1,sd)
hetGP_mean<- apply(hetGPpredMat,1,mean)
hetGP_qMat<- t(as.matrix(apply(hetGPpredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(hetGP_sd,hetGP_mean,hetGP_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGP.RData")


beepr::beep()

NAvals<- apply(hetGPpredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(hetGPpredMat[,i])))
}
sum(NAvals)
which(NAvals>0)

rm(hetGPpredMat)

################################################################################
#homGP no discrepancy

homGPCheappredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/homGPCheap/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  homGPCheappredMat[,i]<- vals
}

homGPCheap_sd<- apply(homGPCheappredMat,1,sd)
homGPCheap_mean<- apply(homGPCheappredMat,1,mean)
homGPCheap_qMat<- t(as.matrix(apply(homGPCheappredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(homGPCheappredMat,homGPCheap_sd,homGPCheap_mean,homGPCheap_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homGPCheap.RData")


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

#hetGPCheap no discrepancy

hetGPCheappredMat<- matrix(NA, nrow= nrow(coords.10m), ncol= nThin)

for(i in 1:nThin){
  run= raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/hetGPCheap/Extent/Run_",i,".asc",sep=""))
  vals<- raster::extract(run,coords.10m)
  hetGPCheappredMat[,i]<- vals
}

hetGPCheap_sd<- apply(hetGPCheappredMat,1,sd)
hetGPCheap_mean<- apply(hetGPCheappredMat,1,mean)
hetGPCheap_qMat<- t(as.matrix(apply(hetGPCheappredMat,1, function(x) quantile(x, probs = c(.05,0.25,.5,.75,.95)) )))


save(hetGPCheappredMat,hetGPCheap_sd,hetGPCheap_mean,hetGPCheap_qMat,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetGPCheap.RData")


beepr::beep()

NAvals<- apply(hetGPCheappredMat,2,function(x) which(is.na(x)))

NAvals<- rep(NA, nThin)
for(i in 1:nThin){
  NAvals[i]<- length(which(is.na(hetGPCheappredMat[,i])))
}
sum(NAvals)
which(NAvals>0)

rm(hetGPCheappredMat)

################################################################################

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/hetML.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/LISFLOOD/homML.RData")

################################################################################