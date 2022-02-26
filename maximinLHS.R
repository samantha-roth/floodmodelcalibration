###Maximin Latin Hypercube Sampling
library(lhs); library(extraDistr)

################################################################################

#Four dimensional Latin Hypercube sample with different priors 
#for river width error and riverbed elevation error

nE<- 200
nCnotE<- 600
nPars= 4

##100 samples, 4 variables
set.seed(22)
lh.E<- lhs::maximinLHS(nE,nPars)
#A 100 by 2 Latin Hypercube Sample matrix with values uniformly distributed on [0,1]
#transform into samples from the right uniform distributions
samp.E<- matrix(0, nE, nPars)
samp.E[,1]<- qunif(lh.E[,1], min = 0.02, max = 0.1)
samp.E[,2]<- qunif(lh.E[,2], min = 0.02, max = 0.4)
samp.E[,3]<- qunif(lh.E[,3], min = 0.95, max = 1.05)
samp.E[,4]<- qunif(lh.E[,4], min = -5, max= 5) 

##300 samples, 4 variables
set.seed(23)
lh.CnotE<- lhs::maximinLHS(nCnotE,nPars)
#A 300 by 2 Latin Hypercube Sample matrix with values uniformly distributed on [0,1]

#transform into samples from the right uniform distributions
samp.CnotE<- matrix(0, nCnotE, nPars)
samp.CnotE[,1]<- qunif(lh.CnotE[,1], min = 0.02, max = 0.1)
samp.CnotE[,2]<- qunif(lh.CnotE[,2], min = 0.02, max = 0.4)
samp.CnotE[,3]<- qunif(lh.CnotE[,3], min = 0.95, max = 1.05) #Uniform (0.95,1.05)
samp.CnotE[,4]<- qunif(lh.CnotE[,4], min = -5, max = 5) #Uniform (-5,5)


samp.C<- rbind(samp.E,samp.CnotE)

colnames(samp.E)<- c("ch","fp","rwe","ree")
colnames(samp.C)<- c("ch","fp","rwe","ree")
colnames(samp.CnotE)<- c("ch","fp","rwe","ree")

save(samp.E,samp.C,samp.CnotE,file=".../FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_allU.RData")

################################################################################
#check samples
load(".../FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_allU.RData")
apply(samp.C,2,summary)

load(".../FloodingModelCalibrationProject/parameterSamples/4param_lhs_samples_tN.RData")
apply(samp.C,2,summary)

