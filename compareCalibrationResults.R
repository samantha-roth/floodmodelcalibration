
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

#identical(res_homGP[2:1e5,],res_homGP[((1e5)+2):(2e5),])
#plot(res_homGP[2:1e5,1],res_homGP[((1e5)+2):(2e5),1])

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML_RW.RData")
#res_homML_RW<- res
#############################################################################################################################

##load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetSML_D.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetSML_D2.RData")
#res_hetMLD<- res

##load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML_D.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML_D2.RData")
#res_homMLD<- res

##load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP_D.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP_D2.RData")
#res_hetGPD<- res

##load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_D.RData")
#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_D2.RData")
#res_homGPD<- res

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP_cheap_D2.RData")
#res_hetGPCheapD<- res

#load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap_D2.RData")
#res_homGPCheapD<- res

#rm(res)
################################################################################

library(batchmeans)

nPars= 5

for(i in 1:nPars){
  print(i)
  #print(paste("HetGP ESS is ",ess(res_hetGP[,i]),sep=""))
  print(paste("HomGP ESS is ",ess(res_homGP[,i]),sep=""))
  print(paste("HomGP800 ESS is ",ess(res_homGP800[,i]),sep=""))
  print(paste("HomML ESS is ",ess(res_homML[,i]),sep=""))
  #print(paste("HetML ESS is ",ess(res_hetML[,i]),sep=""))
  #print(paste("HomGPD ESS is ",ess(res_hetGPD[,i]),sep=""))
  #print(paste("HetGPD ESS is ",ess(res_homGPD[,i]),sep=""))
  #print(paste("HomMLD ESS is ",ess(res_homMLD[,i]),sep=""))
  #print(paste("HetMLD ESS is ",ess(res_hetMLD[,i]),sep=""))
  #print(paste("HomML RW ESS is ",ess(res_homML_RW[,i]),sep=""))
  #print(paste("HetGP (Cheap) ESS is ",ess(res_hetGPCheap[,i]),sep=""))
  print(paste("HomGP (Cheap) ESS is ",ess(res_homGPCheap[,i]),sep=""))
}

for(i in 1:nPars){
  #print(i)
  #acf(res_hetGP[,i],main= paste("HetGP Variable ",i,sep=""))
  acf(res_homGP[,i],main= paste("HomGP Variable ",i,sep=""))
  acf(res_homGP800[,i],main= paste("HomGP800 Variable ",i,sep=""))
  acf(res_homML[,i],main= paste("HomML Variable ",i,sep=""))
  #acf(res_hetML[,i],main= paste("HetML Variable ",i,sep=""))
  #acf(res_hetGPD[,i],main= paste("HetGPD Variable ",i,sep=""))
  #acf(res_homGPD[,i],main= paste("HomGPD Variable ",i,sep=""))
  #acf(res_homMLD[,i],main= paste("HomMLD Variable ",i,sep=""))
  #acf(res_hetMLD[,i],main= paste("HetMLD Variable ",i,sep=""))
  #acf(res_homML_RW[,i],main= paste("HomML RW Variable ",i,sep=""))
  #acf(res_hetGPCheap[,i],main= paste("HetGP (Cheap) Variable ",i,sep=""))
  acf(res_homGPCheap[,i],main= paste("HomGP (cheap) Variable ",i,sep=""))
}

for(i in 1:nPars){
  print(i)
  #print(paste("HetGP est is ",bm(res_hetGP[,i])[1]," and the MCSE is ", bm(res_hetGP[,i])[2], sep=""))
  print(paste("HomGP est is ",bm(res_homGP[,i])[1]," and the MCSE is ", bm(res_homGP[,i])[2], sep=""))
  print(paste("HomGP800 est is ",bm(res_homGP800[,i])[1]," and the MCSE is ", bm(res_homGP[,i])[2], sep=""))
  print(paste("HomML est is ",bm(res_homML[,i])[1]," and the MCSE is ", bm(res_homML[,i])[2], sep=""))
  #print(paste("HetML est is ",bm(res_hetML[,i])[1]," and the MCSE is ", bm(res_hetML[,i])[2], sep=""))
  
  #print(paste("HetGPD est is ",bm(res_hetGPD[,i])[1]," and the MCSE is ", bm(res_hetGPD[,i])[2], sep=""))
  #print(paste("HomGPD est is ",bm(res_homGPD[,i])[1]," and the MCSE is ", bm(res_homGPD[,i])[2], sep=""))
  #print(paste("HomMLD est is ",bm(res_homMLD[,i])[1]," and the MCSE is ", bm(res_homMLD[,i])[2], sep=""))
  #print(paste("HetMLD est is ",bm(res_hetMLD[,i])[1]," and the MCSE is ", bm(res_hetMLD[,i])[2], sep=""))
  #print(paste("HomML RW est is ",bm(res_homML_RW[,i])[1]," and the MCSE is ", bm(res_homML_RW[,i])[2], sep=""))
  #print(paste("HetGP (Cheap) est is ",bm(res_hetGPCheap[,i])[1]," and the MCSE is ", bm(res_hetGPCheap[,i])[2], sep=""))
  print(paste("HomGP (Cheap) est is ",bm(res_homGPCheap[,i])[1]," and the MCSE is ", bm(res_homGPCheap[,i])[2], sep=""))
}


################################################################################
#generate prior samples
priorSampleCh<- runif(1e5,min=.02,max=.1)
priorSampleFp<- runif(1e5,min=.02,max=.4)
priorSampleRWE<- runif(1e5,min=.95,max=1.05)
priorSampleREE<- runif(1e5,min=-5,max=5)

#########CHANNEL ROUGHNESS#######################################################

par(mfrow=c(2,2))
plot(density(res_homML[,1]),main="Channel roughness posterior densities",
     xlab="Channel roughness",
     col="blue",xlim=c(.02, .1), cex.lab = 1.5, cex.axis= 1.5)
lines(density(res_homGP[,1]),col="red")
#lines(density(res_homGP800[,1]),col="turquoise")
#lines(density(res_homML[,1]),col="green")
#lines(density(res_hetGP[,1]),col="yellow")

lines(density(res_homGPCheap[,1]),col="magenta")
#lines(density(res_hetGPCheap[,1]),col="turquoise")
lines(density(priorSampleCh),col="black")
abline(v=.0305,col="gray")

legend(.03,11 , legend=c("HomMR", 
                         "HomGP10",
                         #"HomGP10 (800)",
                         #"HetMR Posterior",
                         #"HetGP Posterior",
                         "HomGP50", 
                         #"HetGP (Cheap) Posterior", 
                         "Prior", "Truth"),
       col=c("blue",
             "red", 
             #"turquoise",
             #"green", 
             #"yellow",
             "magenta",
             #"turquoise",
             "black", "gray"), lty=1,cex=1,box.lty = .1)

#legend(.07, 80, legend=c("HomML Posterior", "HomGP Posterior","HetML Posterior",
#                         "HetGP Posterior",  "HomGP (Cheap) Posterior", 
#                         "HetGP (Cheap) Posterior",  "Prior"),
#       col=c("purple","magenta", "blue", "red","green","yellow", "black"), lty=1,cex=0.8)

###################FLOODPLAIN ROUGHNESS#########################################


plot(density(res_homGP[,2]),main="Floodplain roughness posterior densities",
     xlab="Floodplain roughness",
     col="red",xlim=c(.02, .4), cex.lab = 1.5, cex.axis= 1.5)
lines(density(res_homML[,2]),col="blue")
#lines(density(res_homGP800[,2]),col="turquoise")
#lines(density(res_homML[,2]),col="green")
#lines(density(res_hetGP[,2]),col="yellow")

lines(density(res_homGPCheap[,2]),col="magenta")
#lines(density(res_hetGPCheap[,2]),col="turquoise")
lines(density(priorSampleFp),col="black")
abline(v=.045,col="gray")

legend(.2, 2,legend=c("HomMR", 
                        "HomGP10",
                        #"HomGP10 (800)",
                        #"HomMR Posterior",
                        #"HetGP Posterior",
                        "HomGP50", 
                        #"HetGP (Cheap) Posterior", 
                        "Prior", "Truth"),
       col=c("blue",
             "red", 
             #"turquoise",
             #"green", 
             #"yellow",
             "magenta",
             #"turquoise",
             "black","gray"), lty=1,cex=1,box.lty = 0.1)


#legend(.2, 1.5, legend=c("HomML Posterior", "HomGP Posterior","HetML Posterior",
#                         "HetGP Posterior", "HomGP (Cheap) Posterior", 
#                         "HetGP (Cheap) Posterior", "Prior"),
#       col=c("purple","magenta", "blue", "red","green","yellow","black"), lty=1,cex=0.8)

######################RIVER WIDTH ERROR#########################################


plot(density(res_homGP[,3]),main="River width error posterior densities",
     xlab="River width error",
     col="red",xlim=c(.95, 1.05), cex.lab = 1.5, cex.axis= 1.5)
lines(density(res_homML[,3]),col="blue")
#lines(density(res_homGP800[,3]),col="turquoise")
#lines(density(res_homML[,3]),col="green")
#lines(density(res_hetGP[,3]),col="yellow")

lines(density(res_homGPCheap[,3]),col="magenta")
#lines(density(res_hetGPCheap[,3]),col="turquoise")
lines(density(priorSampleRWE),col="black")
abline(v=1,col="gray")

legend(.99,8 , legend=c("HomMR", 
                      "HomGP10",
                      #"HomGP10 (800)",
                      #"HetMR Posterior",
                      #"HetGP Posterior",
                      "HomGP50", 
                      #"HetGP (Cheap) Posterior", 
                      "Prior","Truth"),
       col=c("blue",
             "red", 
             #"turquoise",
             #"green", 
             #"yellow",
             "magenta",
             #"turquoise",
             "black","gray"), lty=1,cex=1,box.lty = 0.1)


#legend(1, 6, legend=c("HomML Posterior", "HomGP Posterior",
#                      "HetML Posterior","HetGP Posterior", 
#                      "HomGP (Cheap) Posterior", "HetGP (Cheap) Posterior", "Prior"),
#       col=c("purple","magenta", "blue", "red","green","yellow", "black"), lty=1,cex=0.8)



########################RIVERBED ELEVATION ERROR################################
plot(density(res_homGP[,4]),main="Riverbed elevation error posterior densities",
     xlab="Riverbed elevation error",
     col="red",xlim=c(-5,5), cex.lab = 1.5, cex.axis= 1.5)
lines(density(res_homML[,4]),col="blue")
lines(density(res_homGP800[,4]),col="turquoise")
#lines(density(res_homML[,4]),col="green")
#lines(density(res_hetGP[,4]),col="yellow")

lines(density(res_homGPCheap[,4]),col="magenta")
#lines(density(res_hetGPCheap[,4]),col="turquoise")
lines(density(priorSampleREE),col="black")
abline(v=0,col="gray")

legend(-4,.1 , legend=c("HomMR", 
                       "HomGP10",
                       "HomGP10 (800)",
                       #"HomMR Posterior",
                       #"HetGP Posterior",
                       "HomGP50",
                       #"HetGP (Cheap) Posterior", 
                       "Prior","Truth"),
       col=c("blue",
             "red", 
             "turquoise",
             #"green", 
             #"yellow",
             "magenta",
             #"turquoise",
             "black","gray"), lty=1,cex=1,box.lty = 0.1)


#legend(2,.5 , legend=c("HomML Posterior", "HomGP Posterior","HetML Posterior","HetGP Posterior",
#                       "HomGP (Cheap) Posterior", "HetGP (Cheap) Posterior",
#                        "Prior"),
#       col=c("purple","magenta", "blue", "red", "green","yellow", "black"), lty=1,cex=0.8)


#check convergence by looking at the distributions of each parameter for just the first half of the Markov chain then the whole chain


#homML - EH
for(i in 1: ncol(res_homML)){
  plot(density(res_homML[,i]),main=paste("HomML variable ",i,sep=""),col="black")
  lines(density(res_homML[1:round(nrow(res_homML)/2),i]),col="red")
}

#homGP- okay but could go more
for(i in 1: ncol(res_homGP)){
  plot(density(res_homGP[,i]),main=paste("HomGP variable ",i,sep=""),col="black")
  lines(density(res_homGP[1:round(nrow(res_homGP)/2),i]),col="red")
}

#homGP800
for(i in 1: ncol(res_homGP800)){
  plot(density(res_homGP800[,i]),main=paste("HomGP800 variable ",i,sep=""),col="black")
  lines(density(res_homGP800[1:round(nrow(res_homGP)/2),i]),col="red")
}


#homGP Cheap - Fine
for(i in 1: ncol(res_homGPCheap)){
  plot(density(res_homGPCheap[,i]),main=paste("HomGP (Cheap) variable ",i,sep=""),col="black")
  lines(density(res_homGPCheap[1:round(nrow(res_homGPCheap)/2),i]),col="red")
}

#hetML
for(i in 1: ncol(res_hetML)){
  plot(density(res_hetML[,i]),main=paste("HetML variable ",i,sep=""),col="black")
  lines(density(res_hetML[1:round(nrow(res_hetML)/2),i]),col="red")
  #lines(density(res_hetML[1:round(nrow(res_hetML)/5),i]),col="blue")
  #lines(density(res_hetML[1:round(nrow(res_hetML)/10),i]),col="green")
}

#hetGP - NOT FINE
for(i in 1: ncol(res_hetGP)){
  plot(density(res_hetGP[,i]),main=paste("HetGP variable ",i,sep=""),col="black")
  lines(density(res_hetGP[1:round(nrow(res_hetGP)/2),i]),col="red")
}

#hetGP Cheap - Fine
for(i in 1: ncol(res_hetGPCheap)){
  plot(density(res_hetGPCheap[,i]),main=paste("HetGP (Cheap) variable ",i,sep=""),col="black")
  lines(density(res_hetGPCheap[1:round(nrow(res_hetGPCheap)/2),i]),col="red")
}

################################################################################

#hetMLD
for(i in 1: ncol(res_hetMLD)){
  plot(density(res_hetMLD[,i]),main=paste("HetMLD variable ",i,sep=""),col="black")
  lines(density(res_hetMLD[1:round(nrow(res_hetMLD)/2),i]),col="red")
  #ines(density(res_hetMLD[1:round(nrow(res_hetMLD)/5),i]),col="blue")
  #lines(density(res_hetMLD[1:round(nrow(res_hetMLD)/10),i]),col="green")
}

#homMLD
for(i in 1: ncol(res_homMLD)){
  plot(density(res_homMLD[,i]),main=paste("HomMLD variable ",i,sep=""),col="black")
  lines(density(res_homMLD[1:round(nrow(res_homMLD)/2),i]),col="red")
}

#hetGPD
for(i in 1: ncol(res_hetGPD)){
  plot(density(res_hetGPD[,i]),main=paste("HetGPD variable ",i,sep=""),col="black")
  lines(density(res_hetGPD[1:round(nrow(res_hetGPD)/2),i]),col="red")
}

#homGPD
for(i in 1: ncol(res_homGPD)){
  plot(density(res_homGPD[,i]),main=paste("HomGPD variable ",i,sep=""),col="black")
  lines(density(res_homGPD[1:round(nrow(res_homGPD)/2),i]),col="red")
}

#hetGPD Cheap
for(i in 1: ncol(res_hetGPCheapD)){
  plot(density(res_hetGPCheapD[,i]),main=paste("HetGPD (Cheap) variable ",i,sep=""),col="black")
  lines(density(res_hetGPCheapD[1:round(nrow(res_hetGPCheapD)/2),i]),col="red")
}

#homGPD Cheap
for(i in 1: ncol(res_homGPCheapD)){
  plot(density(res_homGPCheapD[,i]),main=paste("HomGPD (Cheap) variable ",i,sep=""),col="black")
  lines(density(res_homGPCheapD[1:round(nrow(res_homGPCheapD)/2),i]),col="red")
}




