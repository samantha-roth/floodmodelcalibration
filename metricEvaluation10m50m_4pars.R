#Process data 
rm(list=ls())
library(boot); library(raster)
library(ggplot2); library(viridis)

nRuns10m<- 200
res.e<-10
res.c<-50
EperC<- (res.c/res.e)^2
#concatenate parameter values
parVals<- matrix(NA,nrow=200,ncol=5)
colnames(parVals)<- c("run","n_ch","n_fp","rwe","ree")

for(i in 1:nRuns10m){
  pars<- read.csv(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/Run_",i,".csv",sep=""))
  parVals[i,]<- as.numeric(pars[,-1])
}
apply(parVals,2,summary)
save(parVals,file="C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")

load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/allParVals.RData")
parVals10m<- as.data.frame(parVals)
#load pseudo observations (run 5)

RunTrue.10m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/Extent/RunTrue_1.asc")
#set crs
crs(RunTrue.10m)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m <- extract(RunTrue.10m,coords.10m)


EuclideanDists10<- rep(0, nRuns10m)
for(i in 1:nRuns10m){
  run<- raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords.10m)
  EuclideanDists10[i]<- sqrt(sum((vals-truevals.10m)^2))
}


#calculate fit and correctness

mGridArea<- res.e^2
rGridArea<- res.e^2

Ar<- length(which(truevals.10m>0))*rGridArea
Fvals10<- rep(0,nRuns10m)
Cvals10<- rep(0,nRuns10m)

for(i in 1:nRuns10m){
  run<- raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs10m/prior1/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords.10m)
  Am<- length(which(vals>0))*mGridArea
  Arm<- length(intersect(which(truevals.10m>0),which(vals>0)))*rGridArea
  Fvals10[i]<- Arm/(Am + Ar - Arm)
  Cvals10[i]<- Arm/Ar
}

#create dataframe of metric values for different parameter values used in different runs
metrics.10m<- cbind(1:length(Fvals10),parVals10m$n_ch,parVals10m$n_fp,parVals10m$rwe,parVals10m$ree,EuclideanDists10,Fvals10,Cvals10)
metrics.10m<- as.data.frame(metrics.10m)
colnames(metrics.10m)<- c("Run","n_ch","n_fp","rwe","ree","EuclideanDists","Fvals","Cvals")

save(metrics.10m,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")

################################################################################
################################################################################

#do this for 50m resolution
nRuns50m<- 800

#concatenate parameter values
parVals<- matrix(NA,nrow=nRuns50m,ncol=5)
colnames(parVals)<- c("run","n_ch","n_fp","rwe","ree")

for(i in 1:nRuns50m){
  pars<- read.csv(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/Run_",i,".csv",sep=""))
  parVals[i,]<- as.numeric(pars[,-1])
}

save(parVals,file="C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")

load("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/allParVals.RData")
parVals50m<- as.data.frame(parVals)
apply(parVals50m,2,summary)

#load pseudo observations (run 5)
RunTrue.50m= raster("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/Extent/RunTrue_1.asc")
crs(RunTrue.50m)<-"+proj=utm +zone=18 +datum=WGS84  +units=m"
coords.50m <- xyFromCell(RunTrue.50m,1:ncell(RunTrue.50m))
truevals.50m<- extract(RunTrue.50m,coords.50m)

EuclideanDists50<- rep(0, nRuns50m)
for(i in 1:nRuns50m){
  run<- raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords.50m)
  EuclideanDists50[i]<- sqrt(sum((vals-truevals.50m)^2))
}


#calculate fit and correctness

mGridArea<- res.c^2
rGridArea<- res.c^2

Ar<- length(which(truevals.50m>0))*rGridArea
Fvals50<- rep(0,nRuns50m)
Cvals50<- rep(0,nRuns50m)

for(i in 1:nRuns50m){
  run<- raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords.50m)
  Am<- length(which(vals>0))*mGridArea
  Arm<- length(intersect(which(truevals.50m>0),which(vals>0)))*rGridArea
  Fvals50[i]<- Arm/(Am + Ar - Arm)
  Cvals50[i]<- Arm/Ar
}

plot(parVals50m$n_ch,Fvals50,xlab= "Channel roughness", ylab= "F",main="Fit of prediction to observation (50m)")
plot(parVals50m$n_ch,Cvals50,xlab= "Channel roughness", ylab= "C",main="Correctness of prediction to observation (50m)")

#create dataframe of metric values for different parameter values used in different runs
metrics.50m<- cbind(1:nRuns50m,parVals50m$n_ch,parVals50m$n_fp,parVals50m$rwe,parVals50m$ree,EuclideanDists50,Fvals50,Cvals50)
metrics.50m<- as.data.frame(metrics.50m)
colnames(metrics.50m)<- c("Run","n_ch","n_fp","rwe","ree","EuclideanDists","Fvals","Cvals")

save(metrics.50m,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics50m.RData")


################################################################################
################################################################################

#what about predicting at the 10m resolution using 50m resolution model output?
plot(RunTrue.50m)
plot(RunTrue.10m)

identical(parVals10m,parVals50m[1:nRuns10m,])

coords.50m <- xyFromCell(RunTrue.50m,1:ncell(RunTrue.50m))
truevals.50m <- extract(RunTrue.50m,coords.50m)

coords.10m <- xyFromCell(RunTrue.10m,1:ncell(RunTrue.10m))
truevals.10m<- extract(RunTrue.10m,coords.10m)

y10m<- unique(coords.10m[,2])
x10m<- unique(coords.10m[,1])
y50m<- unique(coords.50m[,2])
x50m<- unique(coords.50m[,1])

nx10m<- length(x10m)
ny10m<- length(y10m)
nx50m<- length(x50m)
ny50m<- length(y50m)

summary(y10m)
summary(y50m)
summary(x10m)
summary(x50m)


################################################################################
coords.10m<- as.data.frame(coords.10m)
x.u.10<- unique(coords.10m$x)

coords.50m<- as.data.frame(coords.50m)
x.u.50<- unique(coords.50m$x)

coords.10m.round<- apply(coords.10m,2,round)
coords.10m.round<- as.data.frame(coords.10m.round)

min.x.50m<- min(ceiling(x50m))
max.x.50m<- max(ceiling(x50m))
min.y.50m<- min(ceiling(y50m))
max.y.50m<- max(ceiling(y50m))

min.x.10m<- min(ceiling(x10m))
max.x.10m<- max(ceiling(x10m))
min.y.10m<- min(ceiling(y10m))
max.y.10m<- max(ceiling(y10m))

getCorrespInds<- function(xround,yround,indxy){
  logical<- c(xround==min.x.50m,xround==max.x.50m,yround==min.y.50m,yround==max.y.50m)
  
  if(identical(logical,c(FALSE,FALSE,FALSE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2),
                           (indxy-2):(indxy+2),
                           (indxy+nx10m-2):(indxy+nx10m+2),
                           (indxy+2*nx10m-2):(indxy+2*nx10m+2))
  }
  
  if(identical(logical,c(TRUE,FALSE,FALSE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2),
                           (indxy-2):(indxy+2),
                           (indxy+nx10m-2):(indxy+nx10m+2),
                           (indxy+2*nx10m-2):(indxy+2*nx10m+2))
  }
  
  if(identical(logical,c(FALSE,TRUE,FALSE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m),
                           (indxy-nx10m-2):(indxy-nx10m),
                           (indxy-2):(indxy),
                           (indxy+nx10m-2):(indxy+nx10m),
                           (indxy+2*nx10m-2):(indxy+2*nx10m))
  }
  
  if(identical(logical,c(FALSE,FALSE,TRUE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2))
  }
  
  if(identical(logical,c(FALSE,FALSE,FALSE,TRUE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2),
                           (indxy-2):(indxy+2),
                           (indxy+nx10m-2):(indxy+nx10m+2),
                           (indxy+2*nx10m-2):(indxy+2*nx10m+2))
  }
  
  if(identical(logical,c(TRUE,FALSE,TRUE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2))
  }
  
  if(identical(logical,c(TRUE,FALSE,FALSE,TRUE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m+2),
                           (indxy-nx10m-2):(indxy-nx10m+2),
                           (indxy-2):(indxy+2),
                           (indxy+nx10m-2):(indxy+nx10m+2),
                           (indxy+2*nx10m-2):(indxy+2*nx10m+2))
  }
  
  if(identical(logical,c(FALSE,TRUE,TRUE,FALSE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m),
                           (indxy-nx10m-2):(indxy-nx10m))
  }
  
  if(identical(logical,c(FALSE,TRUE,FALSE,TRUE))){
    corresp.coords.inds<-c((indxy-2*nx10m-2):(indxy-2*nx10m),
                           (indxy-nx10m-2):(indxy-nx10m),
                           (indxy-2):(indxy),
                           (indxy+nx10m-2):(indxy+nx10m),
                           (indxy+2*nx10m-2):(indxy+2*nx10m))
  }
  corresp.coords.inds
}

#inds10in50.mat<- matrix(NA, nrow= nx50m*ny50m, ncol= EperC)
#for(i in 1:nrow(coords.50m)){
#  print(i)
#  coords1<- coords.50m[i,]
#  xround<- as.numeric(round(coords1[1]))
#  yround<- as.numeric(round(coords1[2]))
#  indx<- which(coords.10m.round$x==xround)
#  indy<- which(coords.10m.round$y==yround)
#  indxy<- intersect(indx,indy)
  
#  if(yround>min.y.50m){ #indxy contains a point for all y50 vaues above the minimum
#    corresp.coords.inds<- getCorrespInds(xround,yround,indxy)
#  }
#  if(yround==min.y.50m){#indxy is length 0 here
#    yround2<- as.numeric(round(coords1[2]))+res.c
#    indy2<- which(coords.10m.round$y==yround2)
#    indxy2<- intersect(indx,indy2)
#    indxy<- indxy2+ 5*nx10m 
#    #get what the coordinates of the 10m box 
#    #corresponding to the center of the 50m box would be if the 10m box existed
#    corresp.coords.inds<- getCorrespInds(xround,yround,indxy)
#  }
  
#  if(length(corresp.coords.inds)==25){
#    inds10in50.mat[i,]<- corresp.coords.inds
#  }
#  if(length(corresp.coords.inds)<25){
#    inds10in50.mat[i,]<- c(corresp.coords.inds,rep(NA,25-length(corresp.coords.inds)))
#  }

#}
#save(inds10in50.mat, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/inds10in50.mat.RData")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/inds10in50.mat.RData")

#NAvals<- apply(is.na(inds10in50.mat), 2, which)
#numNAs<- 0
#for(j in 1:25){
#  numNAs<- numNAs+ length(NAvals[[j]])
#}

#repeat10mIndices<- rep(NA,nx10m*ny10m)
#which50mbox<- rep(NA,nx10m*ny10m)
#for(i in 1:(nx10m*ny10m)){
#  howmanyappears<- apply(inds10in50.mat,2, function(x) {length(which(x==i))})
#  repeat10mIndices[i]<- sum(howmanyappears)
#  wherevec<- na.omit(as.numeric(apply(inds10in50.mat,2, function(x) {which(x==i)})))
#  if(length(wherevec)>1) print("warning, length > 1")
#  which50mbox[i]<- wherevec
#}

#length(which(repeat10mIndices>1))
#length(which(repeat10mIndices==1)) #all of the indices are covered once as desired

#coords.10m.repeated<- coords.10m[which(repeat10mIndices>1),]
#apply(coords.10m.repeated,2,summary)
#length(which(round(coords.10m[,2])==4516785))

##plot the 50m box number that each 10m box corresponds to in space to make sure everything looks right
#RunTrue.10m.which50mbox<- RunTrue.10m
#values(RunTrue.10m.which50mbox)<- which50mbox
#plot(RunTrue.10m.which50mbox)

#valswhich<- extract(RunTrue.10m.which50mbox,coords.10m)

#original
#yIndsIwant1<- which(coords.10m[,1]>343190)
#yIndsIwant2<- which(coords.10m[,1]<343310)
#xIndsIwant1<- which(coords.10m[,2]>4519050)
#xIndsIwant2<- which(coords.10m[,2]<4519080)
#yIndsIwant<- intersect(yIndsIwant1,yIndsIwant2)
#xIndsIwant<- intersect(xIndsIwant1,xIndsIwant2)


#xIndsIwant<- which(coords.10m[,1]>345204)
#yIndsIwant<- which(coords.10m[,2]>4519770)

#coordsIwantInds<- intersect(yIndsIwant,xIndsIwant)
#coordsIwant<- coords.10m[coordsIwantInds,]

#Run10which50smaller<- rasterFromCells(RunTrue.10m.which50mbox,coordsIwantInds)
#values(Run10which50smaller)<- valswhich[coordsIwantInds]
#plot(Run10which50smaller)
#corresponding 50m x values: 343204.5 343234.5 343264.5 343294.5 

#corresponding 50m y values: 4519065

#corresponding 10m x values: 
#343194.5 343204.5 343214.5 343224.5 343234.5 343244.5 343254.5 343264.5 343274.5 343284.5 343294.5 343504.5

#corresponding 10m y values: 4519075 4519065 4519055

################################################################################

##Compare the coords that were repeated to the regular coords
#summary(coords.10m$x)
#summary(coords.10m.repeated$x)
##all of the x values that were repeated were 341185, the lowest x value
#summary(coords.10m$y)

#summary(coords.10m.repeated$y)
################################################################################
#get predictions

preds10from50<- matrix(NA, nrow= nx10m*ny10m, ncol=nRuns50m)
for(i in 1:nRuns50m){
  if(i%%100==0) print(i)
  run<- raster(paste("C:/FloodingModelCalibrationProject/multires/modelRuns/4Pars/runs50m/prior1/Extent/Run_",i,".asc",sep=""))
  vals<- extract(run,coords.50m)
  for(j in 1:length(vals)){
    inds10<- stats::na.omit(inds10in50.mat[j,])
    preds10from50[inds10,i]<- vals[j] 
  }
}
save(preds10from50,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/preds10from50.RData")

whereNA<- apply(preds10from50,2,function(x){which(is.na(x))}) #no NA values
whereNA

#Euclidean<- function(coords2){
#  coords<- rbind(coords1,coords2)
#  distance<- dist(coords)
#  distance
#}

EuclideanDists10from50<- rep(NA, nRuns50m)
for(i in 1:nRuns50m){
  EuclideanDists10from50[i]<- sqrt(sum((as.numeric(preds10from50[,i])-truevals.10m)^2))
}

#calculate fit and correctness

mGridArea<- res.e^2
rGridArea<- res.e^2

Ar<- length(which(truevals.10m>0))*rGridArea
Fvals<- rep(0,nRuns50m)
Cvals<- rep(0,nRuns50m)

for(i in 1:nRuns50m){
  vals<- as.numeric(preds10from50[,i])
  Am<- length(which(vals>0))*mGridArea
  Arm<- length(intersect(which(truevals.10m>0),which(vals>0)))*rGridArea
  Fvals[i]<- Arm/(Am + Ar - Arm)
  Cvals[i]<- Arm/Ar
}

#sortedFvals<- sort(Fvals,decreasing=TRUE)
#USortedFvals<- unique(sortedFvals)
#indsDecrease<- rep(NA,length(USortedFvals))
#for(i in 1:length(USortedFvals)){indsDecrease[i]<- which(Fvals==USortedFvals[i])[1]}

metrics.10mfrom50m<- cbind(1:length(Fvals),parVals50m$n_ch,parVals50m$n_fp,parVals50m$rwe,parVals50m$ree,EuclideanDists10from50,Fvals,Cvals)
metrics.10mfrom50m<- as.data.frame(metrics.10mfrom50m)
colnames(metrics.10mfrom50m)[1:3]<- c("Run","n_ch","n_fp")
colnames(metrics.10mfrom50m)<- c("Run","n_ch","n_fp","rwe","ree","EuclideanDists","Fvals","Cvals")

save(metrics.10mfrom50m,file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")


################################################################################
################################################################################
############################PLOTS###############################################
################################################################################
################################################################################

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10m.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics50m.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/metrics10mfrom50m.RData")

apply(metrics.10m,2,summary)
apply(metrics.50m,2,summary)
apply(metrics.10mfrom50m,2,summary)

par(mfrow=c(1,3))
plot(metrics.10m$n_ch,metrics.10m$EuclideanDists, xlab="Channel Roughness", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (10m)")
plot(metrics.50m$n_ch,metrics.50m$EuclideanDists, xlab="Channel Roughness", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (50m)")
plot(metrics.10mfrom50m$n_ch,metrics.10mfrom50m$EuclideanDists,xlab="Channel Roughness",ylab="Euclidean Distance",main="Euclidean distance of prediction (50m) from observation (10m)")


par(mfrow=c(1,3))
plot(metrics.10m$n_fp,metrics.10m$EuclideanDists, xlab="Floodplain Roughness", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (10m)")
plot(metrics.50m$n_fp, metrics.50m$EuclideanDists, xlab="Floodplain roughness", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (50m)")
plot(metrics.10mfrom50m$n_fp,metrics.10mfrom50m$EuclideanDists,xlab="Floodplain Roughness",ylab="Euclidean Distance",main="Euclidean distance of prediction (50m) from observation (10m)")


par(mfrow=c(1,3))
plot(metrics.10m$rwe,metrics.10m$EuclideanDists, xlab="River width error", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (10m)")
plot(metrics.50m$rwe, metrics.50m$EuclideanDists, xlab="River width error", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (50m)")
plot(metrics.10mfrom50m$rwe,metrics.10mfrom50m$EuclideanDists, xlab="River width error", ylab="Euclidean Distance",main="Euclidean distance of prediction (50m) from observation (10m)")

par(mfrow=c(1,3))
plot(metrics.10m$ree,metrics.10m$EuclideanDists, xlab="Riverbed elevation error", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (10m)")
plot(metrics.50m$ree, metrics.50m$EuclideanDists, xlab="Riverbed elevation error", ylab="Euclidean Distance",main="Euclidean distance of prediction from observation (50m)")
plot(metrics.10mfrom50m$ree,metrics.10mfrom50m$EuclideanDists, xlab="Riverbed elevation error", ylab="Euclidean Distance",main="Euclidean distance of prediction (50m) from observation (10m)")




library("ggplot2"); library(viridis); library(gridExtra)

#par(mfrow=c(1,3))
ch_fp_10<- ggplot(metrics.10m,aes(x=n_ch,y=n_fp,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")
ch_fp_50<- ggplot(metrics.50m,aes(x=n_ch,y=n_fp,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")
ch_fp_1050<- ggplot(metrics.10mfrom50m,aes(x=n_ch,y=n_fp,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")
grid.arrange(ch_fp_10, ch_fp_50, ch_fp_1050, ncol=3)


#par(mfrow=c(1,3))
ch_rwe_10<- ggplot(metrics.10m,aes(x=n_ch,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("Channel roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
ch_rwe_50<-ggplot(metrics.50m,aes(x=n_ch,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("Channel roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
ch_rwe_1050<-ggplot(metrics.10mfrom50m,aes(x=n_ch,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("Channel roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
grid.arrange(ch_rwe_10, ch_rwe_50, ch_rwe_1050, ncol=3)

#par(mfrow=c(1,3))
ch_ree_10<- ggplot(metrics.10m,aes(x=n_ch,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("Channel roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of channel roughnesses
ch_ree_50<- ggplot(metrics.50m,aes(x=n_ch,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("Channel roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of channel roughnesses
#when channel roughness is too high (low), a lower (higher) riverbed elevation error seems to compensate
ch_ree_1050<- ggplot(metrics.10mfrom50m,aes(x=n_ch,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("Channel roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of channel roughnesses
#when channel roughness is too high (low), a lower (higher) riverbed elevation error seems to compensate
#intersection of dark blue line with 0 ree may be a bit lower than 3
grid.arrange(ch_ree_10, ch_ree_50, ch_ree_1050, ncol=3)

#par(mfrow=c(1,3))
rwe_ree_10<- ggplot(metrics.10m,aes(x=rwe,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("River width error") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too low of river width errors
rwe_ree_50<- ggplot(metrics.50m,aes(x=rwe,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("River width error") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too low of river width errors
rwe_ree_1050<- ggplot(metrics.10mfrom50m,aes(x=rwe,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("River width error") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too low of river width errors
grid.arrange(rwe_ree_10, rwe_ree_50, rwe_ree_1050, ncol=3)


#par(mfrow=c(1,3))
fp_rwe_10<- ggplot(metrics.10m,aes(x=n_fp,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("Floodplain roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
fp_rwe_50<-ggplot(metrics.50m,aes(x=n_fp,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("Floodplain roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
fp_rwe_1050<-ggplot(metrics.10mfrom50m,aes(x=n_fp,y=rwe,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("Floodplain roughness") + ylab("River width error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too low of river width errors and too high of channel roughnesses
grid.arrange(fp_rwe_10, fp_rwe_50, fp_rwe_1050, ncol=3)

#par(mfrow=c(1,3))
fp_ree_10<- ggplot(metrics.10m,aes(x=n_fp,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (10m)") +
  xlab("Floodplain roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of channel roughnesses
fp_ree_50<- ggplot(metrics.50m,aes(x=n_fp,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction from observation (50m)") +
  xlab("Floodplain roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of channel roughnesses
#when channel roughness is too high (low), a lower (higher) riverbed elevation error seems to compensate
fp_ree_1050<- ggplot(metrics.10mfrom50m,aes(x=n_fp,y=ree,col=EuclideanDists)) + geom_point() + ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  xlab("Floodplain roughness") + ylab("Riverbed elevation error") + scale_color_viridis(option = "H")
#Euclidean distances are highest for too high of riverbed elevation errors and too high of fpannel roughnesses
#when fpannel roughness is too high (low), a lower (higher) riverbed elevation error seems to compensate
#intersection of dark blue line with 0 ree may be a bit lower than 3
grid.arrange(fp_ree_10, fp_ree_50, fp_ree_1050, ncol=3)


#ggplot(metrics.10m,aes(x=n_ch,y=n_fp,col=Fvals10)) + geom_point() + ggtitle("Fit of prediction to observation (10m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")

#ggplot(metrics.10m,aes(x=n_ch,y=n_fp,col=Cvals10)) + geom_point() + ggtitle("Correctness of prediction to observation (10m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")

#ggplot(metrics.50m,aes(x=n_ch,y=n_fp,col=Fvals50)) + geom_point() + ggtitle("Fit of prediction to observation (50m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")

#ggplot(metrics.50m,aes(x=n_ch,y=n_fp,col=Cvals50)) + geom_point() + ggtitle("Correctness of prediction to observation (50m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")


#ggplot(metrics.10mfrom50m,aes(x=n_ch,y=n_fp,col=Fvals)) + geom_point() + ggtitle("Fit of prediction (50m) to observation (10m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")

#ggplot(metrics.10mfrom50m,aes(x=n_ch,y=n_fp,col=Cvals)) + geom_point() + ggtitle("Coverage of prediction (50m) to observation (10m)") +
#  xlab("Channel roughness") + ylab("Floodplain roughness") + scale_color_viridis(option = "H")



plot(metrics.10mfrom50m$n_ch,metrics.10mfrom50m$Fvals,xlab= "Channel roughness", ylab= "F")
plot(metrics.10mfrom50m$n_ch,metrics.10mfrom50m$Cvals,xlab= "Channel roughness", ylab= "C")


plot(metrics.10m$n_ch,metrics.10m$Fvals,xlab= "Channel roughness", ylab= "F",main="Fit of prediction to observation (10m)")
plot(metrics.10m$rwe,metrics.10m$Fvals, xlab="River width error", ylab="Euclidean Distance",main="Fit of prediction to observation (10m)")
plot(metrics.10m$ree,metrics.10m$Fvals, xlab="Riverbed elevation error", ylab="Euclidean Distance",main="Fit of prediction to observation (10m)")

plot(metrics.10m$n_ch,metrics.10m$Cvals,xlab= "Channel roughness", ylab= "C",main="Correctness of prediction to observation (10m)")
plot(metrics.10m$rwe,metrics.10m$Cvals, xlab="River width error", ylab="Euclidean Distance",main="Correctness of prediction to observation (10m)")
plot(metrics.10m$ree,metrics.10m$Cvals, xlab="Riverbed elevation error", ylab="Euclidean Distance",main="Correctness of prediction to observation (10m)")


####################################################################################################################

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_ch_ree_10m.jpeg",width = 800, height = 600)
ggplot(metrics.10m,aes(x=n_ch,y=ree,col=EuclideanDists)) + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        legend.text= element_text(size=24)) +
  xlab("Channel roughness") +
  ylab("Riverbed elevation error") + 
  scale_color_viridis(option = "H")
dev.off()

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_ch_ree_10m50m.jpeg",width = 800, height = 600)
ggplot(metrics.10mfrom50m,aes(x=n_ch,y=ree,col=EuclideanDists)) + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        legend.text= element_text(size=24)) +
  xlab("Channel roughness") + 
  ylab("Riverbed elevation error") + 
  scale_color_viridis(option = "H")
dev.off()


#plot(metrics.10m$n_fp,metrics.10m$EuclideanDists, xlab="Floodplain roughness", ylab="Euclidean distance",main="Euclidean distance of prediction from observation (10m)")

jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_fp_10m.jpeg",width = 800, height = 600)
ggplot(metrics.10m,aes(x=n_fp,y=EuclideanDists),col="black") + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (10m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        axis.title = element_text(size=24)) +
  xlab("Floodplain roughness") +
  ylab("Euclidean distance")
dev.off()

#plot(metrics.10m$rwe,metrics.10m$EuclideanDists, xlab="River width error", ylab="Euclidean distance",main="Euclidean distance of prediction from observation (10m)")


jpeg(filename = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_rwe_10m.jpeg",width = 800, height = 600)
ggplot(metrics.10m,aes(x=rwe,y=EuclideanDists),col="black") + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (10m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        axis.title = element_text(size=24)) +
  xlab("River width error") +
  ylab("Euclidean distance")
dev.off()


################################################################################

#saving as pdfs

pdf(file= "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_ch_ree_10m.pdf",width = 800, height = 700)
ggplot(metrics.10m,aes(x=n_ch,y=ree,col=EuclideanDists)) + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        legend.text= element_text(size=24)) +
  xlab("Channel roughness") +
  ylab("Riverbed elevation error") + 
  scale_color_viridis(option = "H")
dev.off()

pdf(file = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_ch_ree_10m50m.pdf",width = 800, height = 700)
ggplot(metrics.10mfrom50m,aes(x=n_ch,y=ree,col=EuclideanDists)) + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (50m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.title = element_text(size=24),
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        legend.text= element_text(size=24)) +
  xlab("Channel roughness") + 
  ylab("Riverbed elevation error") + 
  scale_color_viridis(option = "H")
dev.off()


#plot(metrics.10m$n_fp,metrics.10m$EuclideanDists, xlab="Floodplain roughness", ylab="Euclidean distance",main="Euclidean distance of prediction from observation (10m)")

pdf(file = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_fp_10m.pdf",width = 800, height = 700)
ggplot(metrics.10m,aes(x=n_fp,y=EuclideanDists),col="black") + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (10m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        axis.title = element_text(size=24)) +
  xlab("Floodplain roughness") +
  ylab("Euclidean distance")
dev.off()

#plot(metrics.10m$rwe,metrics.10m$EuclideanDists, xlab="River width error", ylab="Euclidean distance",main="Euclidean distance of prediction from observation (10m)")


pdf(file = "C:/FloodingModelCalibrationProject/multires/plots/4Pars/10mfrom50m/prior1/ED_rwe_10m.pdf",width = 800, height = 700)
ggplot(metrics.10m,aes(x=rwe,y=EuclideanDists),col="black") + 
  geom_point() + 
  ggtitle("Euclidean distance of prediction (10m) from observation (10m)") +
  theme(plot.title = element_text(size=24), 
        axis.text = element_text(size = 20),
        legend.title= element_text(size=0),
        axis.title = element_text(size=24)) +
  xlab("River width error") +
  ylab("Euclidean distance")
dev.off()