
rm(list=ls())
#thinning MCMC samples

#100000 samples

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP.RData")
res_homGP<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homSML.RData")
res_homML<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetSML.RData")
res_hetML<- res[-c(1,(1e5)+1,(2e5)+1 ),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP.RData")
res_hetGP<- res[-c(1,(1e5)+1,(2e5)+1 ),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_hetGP_cheap.RData")
res_hetGPCheap<- res[-c(1,(1e5)+1),]

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_homGP_cheap.RData")
res_homGPCheap<- res[-c(1,(1e5)+1),]

rm(res)
########################################################################################################################

thin_inds<- seq(600, 3e5, by= 600)
#set.seed(52)
#thin_inds<- sample(1:3e5,500,replace=FALSE)

res_homGP_thin<- res_homGP[thin_inds,]
res_hetGP_thin<- res_hetGP[thin_inds,]

res_homGPCheap_thin<- res_homGPCheap[thin_inds,]
res_hetGPCheap_thin<- res_hetGPCheap[thin_inds,]

res_homML_thin<- res_homML[thin_inds,]
res_hetML_thin<- res_hetML[thin_inds,]

save(res_homGPCheap_thin, res_hetGPCheap_thin, res_homGP_thin, res_hetGP_thin, res_homML_thin,res_hetML_thin, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_thin500.RData")

load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_resD_thin500.RData")
#see if thinned densities look like originals

nPars= 5

#hetML
for(i in 1:(nPars-1)){
  plot(density(res_hetML[,i]),main=paste("HetML variable ",i,sep=""),col="black")
  lines(density(res_hetML_thin[,i]),col="red")
  
}

#homML
for(i in 1:(nPars-1)){
  plot(density(res_homML[,i]),main=paste("HomML variable ",i,sep=""),col="black")
  lines(density(res_homML_thin[,i]),col="red")
}

#hetGP
for(i in 1:(nPars-1)){
  plot(density(res_hetGP_thin[,i]),main=paste("HetGP variable ",i,sep=""),col="red")
  lines(density(res_hetGP[,i]),col="black")
}

#homGP
for(i in 1:(nPars-1)){
  plot(density(res_homGP_thin[,i]),main=paste("HomGP variable ",i,sep=""),col="red")
  lines(density(res_homGP[,i]),col="black")
}

#hetGP
for(i in 1:(nPars-1)){
  plot(density(res_hetGPCheap_thin[,i]),main=paste("HetGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_hetGPCheap[,i]),col="black")
}

#homGP
for(i in 1:(nPars-1)){
  plot(density(res_homGPCheap_thin[,i]),main=paste("HomGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_homGPCheap[,i]),col="black")
}

################################################################################

  #generate prior samples
  priorSampleCh<- runif(1e6,min=.02,max=.1)
  priorSampleFp<- runif(1e6,min=.02,max=.4)
  priorSampleRWE<- runif(1e6,min=.95,max=1.05)
  priorSampleREE<- runif(1e6,min=-5,max=5)
  
  #########CHANNEL ROUGHNESS#######################################################
  
  par(mfrow=(c(2,2)))
  
  plot(density(res_hetGPCheap[,1]),main="Channel roughness density",xlab="Channel roughness",col="red",xlim=c(.02, .1))
  lines(density(res_hetML[,1]),col="blue")
  lines(density(res_hetGPCheap_thin[,1]),col="magenta")
  lines(density(res_hetML_thin[,1]),col="turquoise")
  lines(density(priorSampleCh),col="black")
  abline(v=.0305,col="gray")
  
  legend(.03,8, legend=c("HetGP50 Posterior",
                           "HetML Posterior",
                           "HetGP50 Posterior Thinned",
                           "HetML Posterior Thinned", 
                           "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  
  ###################FLOODPLAIN ROUGHNESS#########################################

  
  plot(density(res_hetGPCheap_thin[,2]),main="Floodplain roughness density",xlab="Floodplain roughness",col="magenta",xlim=c(.02, .4))
  lines(density(res_hetML[,2]),col="blue")
  lines(density(res_hetGPCheap[,2]),col="red")
  lines(density(res_hetML_thin[,2]),col="turquoise")
  lines(density(priorSampleFp),col="black")
  abline(v=.045,col="gray")
  
  legend(.15,1.5, legend=c("HetGP50 Posterior",
                         "HetML Posterior",
                         "HetGP50 Posterior Thinned",
                         "HetML Posterior Thinned", 
                         "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  ######################RIVER WIDTH ERROR#########################################
  
  
  
  plot(density(res_hetGPCheap_thin[,3]),main="River width error density",xlab="River width error",col="magenta",xlim=c(.95, 1.05))
  lines(density(res_hetML[,3]),col="blue")
  lines(density(res_hetGPCheap[,3]),col="red")
  lines(density(res_hetML_thin[,3]),col="turquoise")
  lines(density(priorSampleRWE),col="black")
  abline(v=1,col="gray")
  
  legend(.98,8, legend=c("HetGP50 Posterior",
                           "HetML Posterior",
                           "HetGP50 Posterior Thinned",
                           "HetML Posterior Thinned", 
                           "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  ########################RIVERBED ELEVATION ERROR################################
  
  
  plot(density(res_hetGPCheap_thin[,4]),main="Riverbed elevation error density",xlab="Riverbed elevation error",col="magenta",xlim=c(-5, 5))
  lines(density(res_hetML[,4]),col="blue")
  lines(density(res_hetGPCheap[,4]),col="red")
  lines(density(res_hetML_thin[,4]),col="turquoise")
  lines(density(priorSampleREE),col="black")
  abline(v=0,col="gray")
  
  legend(-4,.07, legend=c("HetGP50 Posterior",
                         "HetML Posterior",
                         "HetGP50 Posterior Thinned",
                         "HetML Posterior Thinned", 
                         "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
################################################################################
  

for(i in 1: (ncol(res_hetML)-1)){
  print(i )
  print("HetML thin")
  print(summary(res_hetML_thin[,i]))
  print("HetML")
  print(summary(res_hetML[,i]))
  print("HetGPCheap thin")
  print(summary(res_hetGPCheap_thin[,i]))
  print("HetGPCheap")
  print(summary(res_hetGPCheap[,i]))
}


  
  
################################################################################
  #Hom Plots
  
  #########CHANNEL ROUGHNESS#######################################################
  
  par(mfrow=(c(2,2)))
  
  plot(density(res_homGPCheap[,1]),main="Channel roughness density",xlab="Channel roughness",col="red",xlim=c(.02, .1))
  lines(density(res_homML[,1]),col="blue")
  lines(density(res_homGPCheap_thin[,1]),col="magenta")
  lines(density(res_homML_thin[,1]),col="turquoise")
  lines(density(priorSampleCh),col="black")
  abline(v=.0305,col="gray")
  
  legend(.03,8, legend=c("homGP50 Posterior",
                         "homML Posterior",
                         "homGP50 Posterior Thinned",
                         "homML Posterior Thinned", 
                         "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  
  ###################FLOODPLAIN ROUGHNESS#########################################
  
  
  plot(density(res_homGPCheap_thin[,2]),main="Floodplain roughness density",xlab="Floodplain roughness",col="magenta",xlim=c(.02, .4))
  lines(density(res_homML[,2]),col="blue")
  lines(density(res_homGPCheap[,2]),col="red")
  lines(density(res_homML_thin[,2]),col="turquoise")
  lines(density(priorSampleFp),col="black")
  abline(v=.045,col="gray")
  
  legend(.15,1.5, legend=c("homGP50 Posterior",
                           "homML Posterior",
                           "homGP50 Posterior Thinned",
                           "homML Posterior Thinned", 
                           "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  ######################RIVER WIDTH ERROR#########################################
  
  
  
  plot(density(res_homGPCheap_thin[,3]),main="River width error density",xlab="River width error",col="magenta",xlim=c(.95, 1.05))
  lines(density(res_homML[,3]),col="blue")
  lines(density(res_homGPCheap[,3]),col="red")
  lines(density(res_homML_thin[,3]),col="turquoise")
  lines(density(priorSampleRWE),col="black")
  abline(v=1,col="gray")
  
  legend(.98,8, legend=c("homGP50 Posterior",
                         "homML Posterior",
                         "homGP50 Posterior Thinned",
                         "homML Posterior Thinned", 
                         "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
  ########################RIVERBED ELEVATION ERROR################################
  
  
  plot(density(res_homGPCheap_thin[,4]),main="Riverbed elevation error density",xlab="Riverbed elevation error",col="magenta",xlim=c(-5, 5))
  lines(density(res_homML[,4]),col="blue")
  lines(density(res_homGPCheap[,4]),col="red")
  lines(density(res_homML_thin[,4]),col="turquoise")
  lines(density(priorSampleREE),col="black")
  abline(v=0,col="gray")
  
  legend(-4,.07, legend=c("homGP50 Posterior",
                          "homML Posterior",
                          "homGP50 Posterior Thinned",
                          "homML Posterior Thinned", 
                          "Prior"),
         col=c("red",
               "blue",
               "magenta",
               "turquoise",
               "black"), lty=1,cex=0.8)
  
################################################################################
#check run 147 for homGPD thin and run 148 for hetGPD thin

#res_homGPD_thin[147,]

#res_hetGPD_thin[148,]

#identical_singleRes<- rep(NA, ncol(res_hetGPD))
#for(i in 1:ncol(res_hetGPD)){
#  identical_singleRes[i]<- identical(res_hetGPD[,i],res_homGPD[,i])
#}
#all the same except the last column (variance estimate)


################################################################################
#try to make the thinned densities match the originals better

#thin_inds2<- seq(201, 100001, by= 200)

set.seed(52)
thin_inds2<- sample(1:3e5,500,replace=FALSE)

res_homGP_thin2<- res_homGP[thin_inds2,]
res_hetGP_thin2<- res_hetGP[thin_inds2,]

res_homGPCheap_thin2<- res_homGPCheap[thin_inds2,]
res_hetGPCheap_thin2<- res_hetGPCheap[thin_inds2,]

res_homML_thin2<- res_homML[thin_inds2,]
res_hetML_thin2<- res_hetML[thin_inds2,]


#hetML
for(i in 1: (ncol(res_hetML)-1)){
  par(mfrow=c(1,2))
  plot(density(res_hetML[,i]),main=paste("HetML variable ",i,sep=""),col="black")
  lines(density(res_hetML_thin[,i]),col="red")
  
  plot(density(res_hetML[,i]),main=paste("Thin2: HetML variable ",i,sep=""),col="black")
  lines(density(res_hetML_thin2[,i]),col="red")
  
}

#homML
for(i in 1: (ncol(res_homML)-1)){
  par(mfrow=c(1,2))
  
  plot(density(res_homML[,i]),main=paste("HomML variable ",i,sep=""),col="black")
  lines(density(res_homML_thin[,i]),col="red")
  
  plot(density(res_homML[,i]),main=paste("Thin2: HomML variable ",i,sep=""),col="black")
  lines(density(res_homML_thin2[,i]),col="red")
}

#hetGP
for(i in 1: (ncol(res_hetGP)-1)){
  par(mfrow=c(1,2))
  plot(density(res_hetGP_thin[,i]),main=paste("HetGP variable ",i,sep=""),col="red")
  lines(density(res_hetGP[,i]),col="black")
  
  plot(density(res_hetGP_thin2[,i]),main=paste("Thin2: HetGP variable ",i,sep=""),col="red")
  lines(density(res_hetGP[,i]),col="black")
}

#homGP
for(i in 1: (ncol(res_homGP)-1)){
  par(mfrow=c(1,2))
  plot(density(res_homGP_thin[,i]),main=paste("HomGP variable ",i,sep=""),col="red")
  lines(density(res_homGP[,i]),col="black")
  
  plot(density(res_homGP_thin2[,i]),main=paste("Thin2: HomGP variable ",i,sep=""),col="red")
  lines(density(res_homGP[,i]),col="black")
}

#hetGPCheap
for(i in 1: (ncol(res_hetGPCheap)-1)){
  par(mfrow=c(1,2))
  plot(density(res_hetGPCheap_thin[,i]),main=paste("HetGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_hetGPCheap[,i]),col="black")
  
  plot(density(res_hetGPCheap_thin2[,i]),main=paste("Thin2: HetGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_hetGPCheap[,i]),col="black")
}

#homGPCheap
for(i in 1: (ncol(res_homGPCheap)-1)){
  par(mfrow=c(1,2))
  plot(density(res_homGPCheap_thin[,i]),main=paste("HomGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_homGPCheap[,i]),col="black")
  
  plot(density(res_homGPCheap_thin2[,i]),main=paste("Thin2: HomGP (Cheap) variable ",i,sep=""),col="red")
  lines(density(res_homGPCheap[,i]),col="black")
}

save(res_homGPCheap_thin2, res_hetGPCheap_thin2, res_homGP_thin2, res_hetGP_thin2, res_homML_thin2,res_hetML_thin2, file="C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_thin500.2.RData")
load("C:/FloodingModelCalibrationProject/multires/outputData/4Pars/prior1/10mfrom50m/calibration/mcmc_res_thin500.2.RData")


#see if thinned densities look like originals for a random sample of indices rather than every 200th


################################################################################

#generate prior samples
priorSampleCh<- runif(1e6,min=.02,max=.1)
priorSampleFp<- runif(1e6,min=.02,max=.4)
priorSampleRWE<- runif(1e6,min=.95,max=1.05)
priorSampleREE<- runif(1e6,min=-5,max=5)

#########CHANNEL ROUGHNESS#######################################################

par(mfrow=(c(2,2)))

plot(density(res_hetGPCheapD[,1]),main="Channel roughness density",xlab="Channel roughness",col="red",xlim=c(.02, .1))
lines(density(res_hetMLD[,1]),col="blue")
lines(density(res_hetGPCheapD_thin2[,1]),col="magenta")
lines(density(res_hetMLD_thin2[,1]),col="turquoise")
lines(density(priorSampleCh),col="black")
abline(v=.0305,col="gray")

legend(.03,8, legend=c("HetGP50 Posterior",
                       "HetML Posterior",
                       "HetGP50 Posterior Thinned",
                       "HetML Posterior Thinned", 
                       "Prior"),
       col=c("red",
             "blue",
             "magenta",
             "turquoise",
             "black"), lty=1,cex=0.8)


###################FLOODPLAIN ROUGHNESS#########################################


plot(density(res_hetGPCheapD_thin2[,2]),main="Floodplain roughness density",xlab="Floodplain roughness",col="magenta",xlim=c(.02, .4))
lines(density(res_hetMLD[,2]),col="blue")
lines(density(res_hetGPCheapD[,2]),col="red")
lines(density(res_hetMLD_thin2[,2]),col="turquoise")
lines(density(priorSampleFp),col="black")
abline(v=.045,col="gray")

legend(.15,1.5, legend=c("HetGP50 Posterior",
                         "HetML Posterior",
                         "HetGP50 Posterior Thinned",
                         "HetML Posterior Thinned", 
                         "Prior"),
       col=c("red",
             "blue",
             "magenta",
             "turquoise",
             "black"), lty=1,cex=0.8)

######################RIVER WIDTH ERROR#########################################



plot(density(res_hetGPCheapD_thin2[,3]),main="River width error density",xlab="River width error",col="magenta",xlim=c(.95, 1.05))
lines(density(res_hetMLD[,3]),col="blue")
lines(density(res_hetGPCheapD[,3]),col="red")
lines(density(res_hetMLD_thin2[,3]),col="turquoise")
lines(density(priorSampleRWE),col="black")
abline(v=1,col="gray")

legend(.98,8, legend=c("HetGP50 Posterior",
                       "HetML Posterior",
                       "HetGP50 Posterior Thinned",
                       "HetML Posterior Thinned", 
                       "Prior"),
       col=c("red",
             "blue",
             "magenta",
             "turquoise",
             "black"), lty=1,cex=0.8)

########################RIVERBED ELEVATION ERROR################################


plot(density(res_hetGPCheapD_thin2[,4]),main="Riverbed elevation error density",xlab="Riverbed elevation error",col="magenta",xlim=c(-5, 5))
lines(density(res_hetMLD[,4]),col="blue")
lines(density(res_hetGPCheapD[,4]),col="red")
lines(density(res_hetMLD_thin2[,4]),col="turquoise")
lines(density(priorSampleREE),col="black")
abline(v=0,col="gray")

legend(-4,.07, legend=c("HetGP50 Posterior",
                        "HetML Posterior",
                        "HetGP50 Posterior Thinned",
                        "HetML Posterior Thinned", 
                        "Prior"),
       col=c("red",
             "blue",
             "magenta",
             "turquoise",
             "black"), lty=1,cex=0.8)

################################################################################


#####################################################################################
#####################################################################################
#Hom Plots